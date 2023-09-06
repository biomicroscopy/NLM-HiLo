function opts = set_opts

%--------------- optical parameters
opts.NA_det = 0.5; % detection numerical aperture
opts.NA_ill = 0.5;  % illumination numerical aperture
opts.illumination_wavelength = 480e-9;   %[m]
opts.detection_wavelength = 521e-9;      %[m]
opts.magnification = 8;         % with multi-z relay

%--------------- camera parameters
opts.camera_gain = 2.18;    %[ADU/e-]
opts.camera_noise = 1.2;    %[e-]
opts.pixel_size = 6.5e-6;   %[m]
opts.binning = 1;

%--------------- NLM parameters
opts.applywavelet = 1;
opts.sigma = 1.5;            % HiLo sectioning parameter
opts.block_size = 3;       % block window size [odd] ~ speckle size
opts.search_size = 15;     % search window width = 2*search_size+1    
opts.denoise_flag1 = 1;    % flag for denoising uniform  image
opts.despeckle_flag2 = 1;  % flag for despeckle std  image
opts.h1 = 0.4;             % NLM filtering parameter for denoising u
opts.h2 = 10;              % NLM filtering parameter for despeckle std(u-s)


