%% pre_process_img
% Take acquisition parameters and perform Hilo routine to extract contrast
% weighted uniform image before low pass
% Input:
% Iu - uniform image (offset subtracted)
% Is - speckle image (offset subtracted)
% opts - parameters for reconstruction:
%        opts.sigma: wavelet filter parameter that controls sectioning ability
%        opts.magnification,opts.NA_ill,opts.NA_det: optical parameters
%        opts.pixel_size,opts.binning: camera parameters
%        opts.camera_gain: [ADU/e-]
%        opts.camera_noise: [e-] (rms)

% output:
% Iest - estimate image = C.*Iu
% opts - opts.sigmaHiLo: low pass filtering parameter
%        opts.eta: controls the scaling of lo when merging lo and hi

% 4/25/23 
%%
function [Iest,opts,C] = pre_process_img(Iu,Is,opts)
padnum = 5;
Iu = padarray(Iu,[padnum padnum 0],0);
Is = padarray(Is,[padnum padnum 0],0);

imsize = size(Iu);
sigma_hilo = 2*opts.sigma;


% normalized frequency coordinate [-1/2, 1/2] defined from image size
[F2,Fu,Fv]  = freq_coord_norm(imsize); % normalized frequency coordinate [-1/2, 1/2]
p = opts.pixel_size*opts.binning/opts.magnification;% effective pixel size
K = sqrt(F2)/p; % frequency coordinate

% define camera OTF
OTF_cam = sinc(Fu).*sinc(Fv);
OTF_cam(OTF_cam<0) = 0;

% illumination/detection bandwidth (cutoff = 2NA/wavelength)
bw_ill = 2*opts.NA_ill/opts.illumination_wavelength;
bw_det = 2*opts.NA_det/opts.detection_wavelength;   
% define detection/illumination OTF
OTF_det = real(2/pi*(acos(K/bw_det) - K/bw_det.*sqrt(1-(K/bw_det).^2)));
OTF_det = OTF_det.*OTF_cam;
OTF_ill = real(2/pi*(acos(K/bw_ill) - K/bw_ill.*sqrt(1-(K/bw_ill).^2)));

% define wavelet_filter W(K)
wavelet_gaussian1 = exp(-(pi)^2*F2*opts.sigma^2);
wavelet_gaussian2 = exp(-(pi)^2*F2*2*(opts.sigma)^2);
wavelet_filter = wavelet_gaussian1 - wavelet_gaussian2;
% figure;imagesc(abs(wavelet_filter));axis image;title('Spectrum of wavelet filter');colorbar
% figure;imagesc(real(ft2(wavelet_filter)));axis image;title('Wavelet filter');colorbar

% Normalize
Iu_local = imgaussfilt(Iu,sigma_hilo);
Is_local = imgaussfilt(Is,sigma_hilo);
Iu_filt = Iu./Iu_local;
Is_filt = Is./Is_local;
Id = Iu_filt - Is_filt; % difference image
Id(Iu_filt==0 | Is_filt==0) = 1e-6;

% determine noise
wavelet_volume = sum((abs(wavelet_filter)).^2,'all')/imsize(1)/imsize(2);
noisevar_u = opts.camera_gain.*Iu+(opts.camera_noise*opts.camera_gain).^2;
noisevar_s = opts.camera_gain.*Is+(opts.camera_noise*opts.camera_gain).^2;
noisevar_u = noisevar_u./Iu_local.^2;
noisevar_s = noisevar_s./Is_local.^2;
noisevar =(noisevar_u+noisevar_s)*wavelet_volume*(1-1/sigma_hilo);

% apply wavelet filtering
Id_filt = real(ift2(wavelet_filter.*ft2(Id))); 
C = Id_filt.^2 - noisevar;  
C = (C>=0).*sqrt(C)-(C<0).*sqrt(-C); 
Iest = C.*Iu;

% crop
Iest = Iest(padnum+1:end-padnum,padnum+1:end-padnum,:);

% Calculating eta: scaling factor on "Lo" 
opts.eta = sqrt(sum(OTF_ill(:))/sum(((wavelet_filter.*OTF_det).^2).*OTF_ill,'all'));
opts.sigmaHiLo = sigma_hilo;


function [F2,Fu,Fv] = freq_coord_norm(imsize)
    % Frequency coordinate = fu
    L = imsize;
    fu_line = -1/2 : 1/L(1) : 1/2 - 1/L(1);
    fv_line = -1/2 : 1/L(2) : 1/2 - 1/L(2);
    [Fu,Fv] = meshgrid(fv_line,fu_line);
    F2 = Fu.^2+Fv.^2; 
end



end