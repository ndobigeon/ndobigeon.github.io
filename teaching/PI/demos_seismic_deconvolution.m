clear all
close all
clc

% 1D signal

p = 8;
N = 2^p;

w = 0.05;
z = rand(1,N)<w;
s2 = 10;
x = (sqrt(s2)*randn(1,N)).*z;
x = x(:);

figure(1);
subplot(3,1,1)
plot(1:N,x);
title('spike train')

pause

% wavelet
n_wav = 32;
[wav,supp_wave] = morlet(-3,3,n_wav);

figure(1);
subplot(3,1,2)
plot(supp_wave,wav)
title('wavelet (convolution filter)')

pause

% convolution
y = conv(wav,x);
y = y(n_wav/2:(end-n_wav/2));
figure(1);
subplot(3,1,3)
plot(y);
title('filtered signal')

pause

% matrix formulation
H =  convmtx(wav',N);
y2 = H*x;
size(y2);
y2 = y2(n_wav/2:(end-n_wav/2));

H_trunc = H(n_wav/2:(end-n_wav/2),:);

figure(2);
imagesc(H_trunc);

figure(1)
hold on
plot(y2,'r*')

pause

% naive deconvolution
x_hat = H_trunc\y;

figure;
subplot(311)
plot(1:N,x)
title('spike train')
subplot(312)
plot(1:N,y)
title('filtered signal')
subplot(313)
plot(1:N,x_hat)
title('deconvolved signal')


pause

% noisy data
SNRdB = 80;
Px = 1/N*norm(x)^2;
sigma2 = 10^(-SNRdB/10)*Px;
b = sqrt(sigma2)*randn(N,1);
y_noisy = y+b;



% naive deconvolution
x_hat_noisy = H_trunc\y_noisy;

figure;
subplot(311)
plot(1:N,x)
title('spike train')
subplot(312)
plot(1:N,y_noisy)
hold on
plot(1:N,y,'g')
hold off
title('(noisy) filtered signal')
pause
subplot(313)
plot(1:N,x_hat_noisy)
title('deconvolved noisy signal')

