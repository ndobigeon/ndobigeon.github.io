clear all
close all
clc

%% Generation of the input signal

% 1D input (unknown) signal
p = 8;
N = 2^p;

w = 0.05;
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

% noisy data
SNRdB  = 80;
Px     = 1/N*norm(x)^2;
sigma2 = 10^(-SNRdB/10)*Px;
b = sqrt(sigma2)*randn(N,1);
y_noisy = y+b;

%% Naive inversion
% do not compute H^(-1)
% use left-division \ since it is more computationally efficient

x_naive = H\y_noisy;

%% Resolution by OMP
% do not compute H^(-1)

thres = sqrt(sigma2*N)*10;
T = N;
[x_omp, Tab_a, Tab_k, Tab_r] = algo_OMP(y,H,thres,T);


%% plots

% evolution along the MP iterations

T = length(Tab_a);
xhat = zeros(N,1);

figure(2)
pause
for t=1:T
    xhat(Tab_k(1:t)) = Tab_a{t};
    yhat = H*xhat;
    figure(2)
    subplot(3,1,1)
    plot(1:N,y,'m');
    hold on
    plot(1:N,y_noisy);
    plot(1:N,yhat,'k');
    hold off
    axis tight
    title(['output signal y(n) with SNR=' num2str(SNRdB) 'dB'])
    legend('noise-free signal','noisy signal', 'reconstructed signal')
    subplot(3,1,2)
    plot(1:N,x);
    hold on
    plot(1:N,xhat);
    hold off
    legend('true input signal','estimated signal')
    title('estimated signal - OMP')
    axis tight
    subplot(3,1,3)
    plot(1:T,[Tab_r(1:t) NaN*zeros(1,T-t)])
    axis([1 T 0 max(Tab_r)])
    title('residual energy')    
    pause(0.01)
end


figure(1);
subplot(4,1,1)
plot(1:N,x);
title('input signal x(n) [spike train]')
axis tight
subplot(4,1,2)
plot(1:N,y,'m');
hold on
plot(1:N,y_noisy);
hold off
legend('noise-free signal','noisy signal')
title(['output signal y(n) with SNR=' num2str(SNRdB) 'dB'])
axis tight
subplot(4,1,3)
plot(1:N,x);
hold on
plot(1:N,x_naive,'k-');
hold off
legend('true input signal','estimated signal')
title('estimated signal - naive inversion')
axis tight
subplot(4,1,4)
plot(1:N,x);
hold on
plot(1:N,x_omp,'k-');
hold off
legend('true input signal','estimated signal')
title('estimated signal - OMP')
axis tight


