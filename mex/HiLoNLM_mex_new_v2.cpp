/*=================================================================
% denoise - denoise an image.
%
%   [M1,Wx,Wy] = perform_nlmeans_vectorized(Ma,H,Ha,Vx,Vy,T,max_dist,do_median,mask,exlude_self);
%	[img_std,img_u_new] = HiLoNLM(img_std,Hu,Hulogu_sum,block_size,search_size,lambda);
%
%	Ma is the image used to perform denoising.
%	H is a high dimensional representation (generaly patch wise) of the image to denoise.
%	Ha is a h.d. representation for Ma.
%	Vx is the x position of center of the search in Ma (same size as H).
%	Vy is the y position of center of the search in Ma (same size as H).
%	T is the variance of the gaussian used to compute the weights.
%	max_dist restricts the search size around position given by Vx,Vy.
%	mask_process restrict the area of filtering (useful for inpainting)
%	exclude_self avoid to take into account the central pixel
%
%	M1 is the denoised image
%	Wx is the new center x position for next seach (center of best fit)
%	Wy is the new center y position for next seach (center of best fit)
%
%   Copyright (c) 2006 Gabriel Peyr?
%	11/27/22 update mexFunction
%   11/28/22 dist_windows
%   12/19/22 debugged C++ update
%   01/11/23 seperate realizations of denoise and despeckle 
*=================================================================*/

#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <string.h>
//#include <time.h>
#include <ppl.h>

using namespace concurrency;

// preprocessing phase, define macros,executed before the source code is compiled
#define accessH(A,a,b,c) A[(a)+m2*(b)+m2*n2*(c)]
#define accessSearchinH(B,a,b,c) B[(a)+block_num1D*(b)+block_num1D*block_num1D*(c)]
#define accessImgstd(A,a,b,c) A[(a)+m1*(b)+m1*n1*(c)]
#define accessOutput(C,a,b,c) C[(a)+m0*(b)+m0*n0*(c)]

#define Hu_(a,b,c) accessH(Hu,a,b,c)
#define Hulogu_sum_(a,b) accessH(Hulogu_sum,a,b,0)
#define w_(a,b) accessSearchinH(w,a,b,0)
#define img_std_(a,b) accessImgstd(img_std,a,b,0)
#define output_std_(a,b,c) accessOutput(output_std,a,b,c)
#define output_u_(a,b,c) accessOutput(output_u,a,b,c)



/* Global variables */
int n1 = -1;				// width  of img_std
int m1 = -1;				// height of img_std
int n2 = -1;				// width  of Hu
int m2 = -1;				// height of Hu
int n0 = -1;				// width  of output
int m0 = -1;				// height of output
int s = -1;					// number of color chanels of M,Ma
int k1 = -1;					// dimensionality of the vectorized patches H,Ha

float* output_std = NULL;	// output
float* output_u = NULL;		// output
float* img_std = NULL;	    // input std image of u-s 
float* Hu = NULL;			// vectorized patches from img_u
float* Hulogu_sum = NULL;	// vectorized patches from u*log(u), summed along 3rd dimension
float* w = NULL;			// weight matrix

int block_size = 3;			// patch size
int search_size = 10;		// search window half size
float lambda1 = 10.0f;		// width of the gaussian for denoise
float lambda2 = 10.0f;		// width of the gaussian for despeckle

int block_num1D = (search_size * 2 + 1) - block_size + 1;
int H_search_size = block_num1D / 2;
int half_block_size = block_size / 2;
bool exclude_self = false;
bool denoise_u = false;
bool despeckle_std = false;

inline void display_message(const char* mess, int v)
{
	char str[128];
	sprintf(str, mess, v);
	mexPrintf(str);
	mexPrintf(",\t");
}
inline void display_messagef(const char* mess, float v)
{
	char str[128];
	sprintf(str, mess, v);
	mexPrintf(str);
	mexPrintf(",\t");
}

//Pointer to an mxArray array and get array dimensions [a, b, c]
inline void get_dimensions(const mxArray* arr, int& a, int& b, int& c)
{
	int nd = mxGetNumberOfDimensions(arr);
	if (nd == 3)
	{
		a = mxGetDimensions(arr)[0];
		b = mxGetDimensions(arr)[1];
		c = mxGetDimensions(arr)[2];
	}
	else if (nd == 2)
	{
		a = mxGetM(arr);
		b = mxGetN(arr);
		c = 1;
	}
	else
		mexErrMsgTxt("only 2D and 3D arrays supported.");
}



/* the distance function used */
inline float weight_glr(float x, float LAMBDA)
{
	// x is the GLR
	return exp(x / (k1 * LAMBDA));
}


//// #define weight_func 
//#define weight_func weight_glr

/*
	compute the general likelihood ratio between two windows.
	(i,j) is the center for the reference patch to denoise (in Hu)
	(i1,j1) is the center for the surrounding patches (in Hu),
	p is the dimensionality of the vectors (#pixels in patch)
	GLR = sum[(x1+x2)(log(x1+x2)-log(2)) - x1log(x2) - x2log(x2)]
	x1 and x2 are pixels in each patch
*/
inline
float dist_windows(int i, int j, int i1, int j1)
{
	float dist_GLR = 0;
	for (int a = 0; a < k1; ++a)
	{
		float d = Hu_(i, j, a) + Hu_(i1, j1, a);
		d = d * (log(d) - 0.6931f);
		dist_GLR += d;
	}
	dist_GLR = dist_GLR - Hulogu_sum_(i, j) - Hulogu_sum_(i1, j1);
	return dist_GLR;
}

/*
	compute weights for each patch in the search window
	reference patch centered at pixel (i,j)
	search window range defined by i_min, i_max, j_min, j_max
	weight matrix w has size: [block_num1D, block_num1D]
*/

//float compute_weights(int i, int j, int i_min, int i_max, int j_min, int j_max)
//{
//	float w_sum = 0;	// sum of all weights
//	for (int i1 = i_min; i1 <= i_max; ++i1)	// loop through search window
//		for (int j1 = j_min; j1 <= j_max; ++j1)
//		{
//			if ((!exclude_self) || (i1 != i) || (j1 != j))
//			{
//				float ww = dist_windows(i, j, i1, j1);
//				ww = weight_func(ww);
//				w_sum += ww;
//				w_(i1 - i_min, j1 - j_min) = ww; // store weights
//			}
//			else
//				w_(i1 - i_min, j1 - j_min) = 0;
//		}
//	return w_sum;
//}


void denoise()
{
	display_message("block size=%d", block_size);
	display_message("half search size=%d", search_size);
	display_messagef("(denoise only) lambda1=%g", lambda1);

	parallel_for(int(0), m0, [&](int i)
		{
			//for (int i = 0; i < m0; ++i) 
			for (int j = 0; j < n0; ++j)
			{
				// center for the search in H
				int ic = i + H_search_size;
				int jc = j + H_search_size;
				// compute search region around center pixel in H
				int i_min = ic - H_search_size;
				int i_max = ic + H_search_size;
				int j_min = jc - H_search_size;
				int j_max = jc + H_search_size;
				/* perform reconstruction of std*/
				for (int a = 0; a < s; ++a)		//channel number s = 1
				{
					output_std_(i, j, a) = 0;
					output_u_(i, j, a) = 0;
					float w_sum1 = 0;	// sum of all weights
					float ww, ww1;
					for (int i1 = i_min; i1 <= i_max; ++i1)	// loop through search window
						for (int j1 = j_min; j1 <= j_max; ++j1)
						{
							// compute weight for each patches in search window
							float ww;
							if ((!exclude_self) || (i1 != i) || (j1 != j))
							{
								ww = dist_windows(ic, jc, i1, j1);
							}
							else {
								ww = 0.0f;
							}
							ww1 = weight_glr(ww, lambda1);
							w_sum1 += ww1;
							output_u_(i, j, a) = output_u_(i, j, a) + ww1 * Hu_(i1, j1, k1 / 2);

						}
						output_u_(i, j, a) = output_u_(i, j, a) / w_sum1;

				}
			}
			});

	//mexPrintf("%d\n", i);
	}

void despeckle()
{
	display_message("block size=%d", block_size);
	display_message("half search size=%d", search_size);
	display_messagef("(despeckle only) lambda2=%g", lambda2);

	parallel_for(int(0), m0, [&](int i)
		{
			//for (int i = 0; i < m0; ++i) 
			for (int j = 0; j < n0; ++j)
			{
				// center for the search in H
				int ic = i + H_search_size;
				int jc = j + H_search_size;
				// compute search region around center pixel in H
				int i_min = ic - H_search_size;
				int i_max = ic + H_search_size;
				int j_min = jc - H_search_size;
				int j_max = jc + H_search_size;
				/* perform reconstruction of std*/
				for (int a = 0; a < s; ++a)		//channel number s = 1
				{
					output_std_(i, j, a) = 0;
					output_u_(i, j, a) = 0;
					float w_sum2 = 0;	// sum of all weights
					float ww, ww2;
					for (int i1 = i_min; i1 <= i_max; ++i1)	// loop through search window
						for (int j1 = j_min; j1 <= j_max; ++j1)
						{
							// compute weight for each patches in search window
							float ww;
							if ((!exclude_self) || (i1 != i) || (j1 != j))
							{
								ww = dist_windows(ic, jc, i1, j1);
							}
							else {
								ww = 0.0f;
							}
							ww2 = weight_glr(ww, lambda2);
							w_sum2 += ww2;
							output_std_(i, j, a) = output_std_(i, j, a) + ww2 * img_std_(i1 + half_block_size, j1 + half_block_size);

						}
					output_std_(i, j, a) = output_std_(i, j, a) / w_sum2;

				}
			}
		});
}


void denoise_n_speckle()
{
	display_message("block size=%d", block_size);
	display_message("half search size=%d", search_size);
	display_messagef("(both enabled ) lambda1=%g", lambda1);
	display_messagef("lambda2=%g", lambda2);


	parallel_for(int(0), m0, [&](int i)
		{
			//for (int i = 0; i < m0; ++i) 
			for (int j = 0; j < n0; ++j)
			{
				// center for the search in H
				int ic = i + H_search_size;
				int jc = j + H_search_size;
				// compute search region around center pixel in H
				int i_min = ic - H_search_size;
				int i_max = ic + H_search_size;
				int j_min = jc - H_search_size;
				int j_max = jc + H_search_size;
				/* perform reconstruction of std*/
				for (int a = 0; a < s; ++a)		//channel number s = 1
				{
					output_std_(i, j, a) = 0;
					output_u_(i, j, a) = 0;
					float w_sum1 = 0, w_sum2 = 0;	// sum of all weights
					float ww, ww1, ww2;
					for (int i1 = i_min; i1 <= i_max; ++i1)	// loop through search window
						for (int j1 = j_min; j1 <= j_max; ++j1)
						{
							// compute weight for each patches in search window
							float ww;
							if ((!exclude_self) || (i1 != i) || (j1 != j))
							{
								ww = dist_windows(ic, jc, i1, j1);
							}
							else {
								ww = 0.0f;
							}
							ww1 = weight_glr(ww, lambda1);
							w_sum1 += ww1;
							output_u_(i, j, a) = output_u_(i, j, a) + ww1 * Hu_(i1, j1, k1 / 2);

							ww2 = weight_glr(ww, lambda2);
							w_sum2 += ww2;
							output_std_(i, j, a) = output_std_(i, j, a) + ww2 * img_std_(i1 + half_block_size, j1 + half_block_size);

						}
					output_u_(i, j, a) = output_u_(i, j, a) / w_sum1;
					output_std_(i, j, a) = output_std_(i, j, a) / w_sum2;

				}
			}
		});

}

		void mexFunction(int nlhs, mxArray * plhs[],
			int nrhs, const mxArray * prhs[])
	{
		if (nrhs < 4)
			mexErrMsgTxt("Minimum 4 input arguments required.");
		if( nlhs > 2 )
		 	mexErrMsgTxt("Maximum 2 output arguments.");

		// -- input 1 : img_std -- (m1 = m - block_size + 1)
		get_dimensions(prhs[0], m1, n1, s);
		img_std = (float*)mxGetData(prhs[0]);

		// -- input 2 : Hu -- (k = block_size^2)
		get_dimensions(prhs[1], m2, n2, k1);
		Hu = (float*)mxGetData(prhs[1]);

		// -- input 3 : Hulogu_sum --
		get_dimensions(prhs[2], m2, n2, s);
		Hulogu_sum = (float*)mxGetData(prhs[2]);

		// -- input 4 : block_size --
		block_size = (int)*mxGetPr(prhs[3]);

		// -- input 5 : search_size --
		if (nrhs >= 5)
			search_size = (int)*mxGetPr(prhs[4]);

		// -- input 6 : denoise_u --
		denoise_u = false;
		if (nrhs >= 6)
			denoise_u = (int)*mxGetPr(prhs[5]) > 0.5;

		// -- input 7 : lambda1 --
		if (nrhs >= 7)
			lambda1 = (float)*mxGetPr(prhs[6]);

		// -- input 8 : despeckle_std --
		despeckle_std = false;
		if (nrhs >= 8)
			despeckle_std = (int)*mxGetPr(prhs[7]) > 0.5;

		// -- input 9 : lambda2 --
		if (nrhs >= 9)
			lambda2 = (float)*mxGetPr(prhs[8]);

		if (nlhs > 10)
			mexErrMsgTxt("Too many input arguments.");

		m0 = m1 - 2 * search_size; // raw output image size
		n0 = n1 - 2 * search_size; // raw output image size

		block_num1D = (search_size * 2 + 1) - block_size + 1;
		H_search_size = block_num1D / 2;
		half_block_size = block_size / 2;

		// -- outpout 1 : output_std -- 
		mwSize dims[3] = { m0,n0,s };
		plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
		output_std = (float*)mxGetData(plhs[0]);

		// -- outpout 2 : image_u_new -- 
		plhs[1] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
		output_u = (float*)mxGetData(plhs[1]);

		// matrix of the weights for each pixel 
		w = (float*)malloc(block_num1D * block_num1D * sizeof(float));
		memset(w, 0, block_num1D * block_num1D * sizeof(float));

		/* Do the actual computations in a subroutine */
		if (denoise_u && !despeckle_std) {
			denoise();
		}
		if (!denoise_u && despeckle_std) {
			despeckle();
		}
		if (denoise_u && despeckle_std) {
			denoise_n_speckle();
		}
		free(w);
	}



	// structure to hold result of qsort
	//struct pixel
	//{
	//	double v;
	//	int i;
	//	int j;
	//};
	//
	//// compare a1 and a2 pointed pixels
	//int pixel_cmp(const void* a1, const void* a2)
	//{
	//	const pixel* c1 = (const pixel*)a1;
	//	const pixel* c2 = (const pixel*)a2;
	//	if (c1->v < c2->v)
	//		return -1;
	//	else if (c1->v > c2->v)
	//		return 1;
	//	return 0;
	//}
		// -- input 12 : exclude_self
	//	exclude_self = false;
	//	if (nrhs >= 12)
	//		exclude_self = *mxGetPr(prhs[11]) > 0.5;
	//	// check

		//inline void display_messagef(const char* mess, float v)
//{
//	char str[128];
//	sprintf(str, mess, v);
//	mexWarnMsgTxt(str);
//}