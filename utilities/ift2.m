function b=ift2(a)

% b=fftshift(fft2(a));

b =  ifftshift(ifft2(ifftshift(a)));