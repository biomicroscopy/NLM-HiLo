%% hilo reconstruction with NLM 
% INPUT:
% - img_est0: estimate image
% - u: uniform image
% -opts: NLM options
% Output:
% - output_est: de-speckled estimate image
% - output_u: denoised uniform image
% 12/21/2022
%%
function [output_est, output_u] = hilo_mex_recon(img_est0, u, opts)
block_size = opts.block_size;
search_size = opts.search_size;

%% Re-organize input matrices for mex function
% pad images
img_u = single(padarray(u, [search_size,search_size],'symmetric','both'));
img_est = single(padarray(img_est0, [search_size,search_size],'symmetric','both'));
[m,n] = size(img_u);

% pre-calculate high dimesion representation of u*log(u) and u [single frame]
Hu = im2col(img_u,[block_size,block_size],'sliding');
Hu = shiftdim(reshape(Hu,size(Hu,1),m-block_size+1,[]),1);% vectorized space
ulogu = img_u.*log(img_u);
% Hulogu_sum = conv2(ulogu,ones(block_size),'valid');
Hulogu_sum = im2col(ulogu,[block_size,block_size],'sliding');
Hulogu_sum = shiftdim(reshape(Hulogu_sum,size(Hulogu_sum,1),m-block_size+1,[]),1);
Hulogu_sum = sum(Hulogu_sum,3);

%% mex
[output_est,output_u] = HiLoNLM_mex_new_v2(img_est,Hu,Hulogu_sum,block_size,search_size,opts.denoise_flag1,opts.h1,opts.despeckle_flag2,opts.h2);

if opts.despeckle_flag2==0
    output_est = img_est0;
end

if opts.denoise_flag1==0
    output_u = u;
end

