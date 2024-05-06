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

% noisy data
SNRdB  = 45;
Px     = 1/N*norm(x)^2;
sigma2 = 10^(-SNRdB/10)*Px;
b = sqrt(sigma2)*randn(N,1);
y_noisy = y+b;


%% LS inversion
% pseudo inverse
x_LS = H\y_noisy;

figure(1);
subplot(3,1,1)
plot(1:N,x);
title('input signal x(n) [spike train]')
axis tight
subplot(3,1,2)
plot(1:N,y,'m');
hold on
plot(1:N,y_noisy);
hold off
legend('noise-free signal','noisy signal')
title(['output signal y(n) with SNR=' num2str(SNRdB) 'dB'])
axis tight
subplot(3,1,3)
plot(1:N,x);
hold on
plot(1:N,x_LS,'k-');
hold off
legend('true input signal','LS solution')
title('estimated signal')
axis tight


%% TSVD inversion

[U,S,V] = svd(full(H));
U=U';
S = diag(S);

figure(2)
semilogy(1:size(S,1),S)
xlabel('eigenvalue index')
ylabel('eigenvalue')

% condition number
cond(H)
max(S)/min(S)

Tab_MSE = zeros(1,N)+NaN;
Tab_fit = zeros(1,N)+NaN;

figure(3)
for K = 1:1:N

    x_TSVD = V(:,1:K)*(diag(1./S(1:K))*U(1:K,:)*y_noisy);
        
    Tab_MSE(K) = norm(x-x_TSVD)^2;
    Tab_fit(K) = norm(y_noisy-H*x_TSVD);
    
    subplot(311)
    plot(1:N,x)
    axis tight
    hold on
    plot(1:N,x_TSVD,'k-')
    hold off
    title(['K = ' num2str(K)])
    subplot(312)
    semilogy(1:N,Tab_MSE)
    xlabel('K')
    ylabel('$||x-\hat{x}_{\mathrm{TSVD}}^{(K)}||_2^2$','interpreter','latex')
    xlim([1 N]);
    subplot(313)
    semilogy(1:N,Tab_fit)
    xlabel('K')
    ylabel('$||y-H\hat{x}_{\mathrm{TSVD}}^{(K)}||_2^2$','interpreter','latex')
    xlim([1 N]);
    drawnow
end


% looking for optimal K
[minerr, K_opt] = min(Tab_MSE);
x_TSVD_opt = V(:,1:K_opt)*(diag(1./S(1:K_opt))*U(1:K_opt,:)*y_noisy);

figure(4)
semilogy((Tab_MSE),'k')
hold on
plot(K_opt,Tab_MSE(K_opt),'ro')
text('string',['K_{opt} = ' num2str(K_opt) ],'units','normalized','position',[0.61 0.08],'color','r')
xlabel('K')
ylabel('$||x-\hat{x}_{\mathrm{TSVD}}^{(K)}||_2^2$','interpreter','latex')

figure(5);
subplot(3,1,1)
plot(1:N,x);
title('input signal x(n) [spike train]')
axis tight
subplot(3,1,2)
plot(1:N,y,'m');
hold on
plot(1:N,y_noisy);
hold off
legend('noise-free signal','noisy signal')
title(['output signal y(n) with SNR=' num2str(SNRdB) 'dB'])
axis tight
subplot(3,1,3)
plot(1:N,x);
hold on
plot(1:N,x_LS,'y-');
plot(1:N,x_TSVD_opt,'k-');
axis([1 N min(x_TSVD_opt)*1.1 max(x_TSVD_opt)*1.1])
hold off
legend('true input signal','LS solution','TSVD solution')
title('estimated signal')



