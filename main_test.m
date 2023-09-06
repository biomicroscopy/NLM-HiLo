addpath('utilities\')
addpath('mex\')

%% Load data 
load('test_data.mat');
offset = 118;
u = double(u) - offset; s = double(s) - offset;
u(u<=0) = 1e-6; s(s<=0) = 1e-6;

figure(1);
subplot(1,2,1);
imagesc(u(:,:,1));axis image;colormap gray;title('Uniform illumination')
subplot(1,2,2);
imagesc(s(:,:,1));axis image;colormap gray;title('Speckle illumination')

%% NLM HiLo 
opts = set_opts;
[img_est0,opts] = pre_process_img(u,s,opts);
output_est = zeros(size(u));
output_u = zeros(size(u));

for z = 1%:size(u,3)  
    [output_est_temp, output_u_temp] = hilo_mex_recon(img_est0(:,:,z), u(:,:,z), opts);
    output_est(:,:,z) = output_est_temp;
    output_u(:,:,z) = output_u_temp;
end
disp('NLM done')

% % combine Hi and Lo [NLM]
hi = output_u-imgaussfilt(output_u,opts.sigmaHiLo);
lo = imgaussfilt(output_est,opts.sigmaHiLo); lo(lo<0) = 1e-6;

eta = estimate_eta(hi,lo);
hilo = hi + eta.*lo; hilo(hilo<0)=0;

%% Basic HiLo
hi0 = u-imgaussfilt(u,opts.sigmaHiLo);
lo0 = imgaussfilt(img_est0,opts.sigmaHiLo); lo(lo<0) = 1e-6;
hilo0 = hi0 + eta.*lo0;  hilo0(hilo0<0)=0;

%% Display
figure(21);
subplot(1,2,1);
imagesc(hilo0(:,:,1));axis image;colormap gray;title('Basic HiLo')
subplot(1,2,2);
imagesc(hilo(:,:,1));axis image;colormap gray;title('NL-means HiLo')



