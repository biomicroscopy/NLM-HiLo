function b=ft2(a)

% b=ifft2(ifftshift(a));

b= fftshift(fft2(fftshift(a)));
