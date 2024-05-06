clear all
close all
clc

%% Generation of the input signal

% 1D input (unknown) signal
p = 8;
N = 2^p;

w = 0.1;
a = rand(1,N)<w;
s2 = 10;
x = (sqrt(s2)*randn(1,N)).*a;
x = x(:);

%% Modeling of the forward model (FIR filter with Morlet IR)

% wavelet
n_wav = 32;
P = n_wav-1;
[wav,supp_wave] = morlet(-3,3,n_wav);

%% Generation of the output signal

% convolution by matrix formulation
Hfull = convmtx(wav',N);
H     = Hfull(n_wav/2:(end-n_wav/2),:);
y     = H*x;

%% Naive inversion
% do not compute H^(-1)
% use left-division \ since it is more computationally efficient

xhat = H\y;

figure(1);
subplot(3,1,1)
plot(1:N,x);
title('input signal x(n) [spike train]')
axis tight
subplot(3,1,2)
plot(1:N,y);
title('output signal y(n)')
axis tight
subplot(3,1,3)
plot(1:N,x);
hold on
plot(1:N,xhat,'k.');
hold off
title('estimated signal')
axis tight
