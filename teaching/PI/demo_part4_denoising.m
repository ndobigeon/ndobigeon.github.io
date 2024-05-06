clear all
close all
clc
isprint = 1;
graphicpath = '../slides/figures/';
    
%% Generation of the unknown signal

% 1D input (unknown) signal
p = 9;
N = 2^p;
vect_t = (1:N)';
fe = 500;
f0 = 1;
f1 = 30;
x = cos(2*pi*f0/fe*vect_t+pi/2) + 1/4*cos(2*pi*f1/fe*vect_t);

%% Generation of the measurements

% noisy data
SNRdB  = 20;
Px     = 1/N*norm(x)^2;
sigma2 = 10^(-SNRdB/10)*Px;
b = sqrt(sigma2)*randn(N,1);
y = x+b;

% spectral analysis - periodogram
Sx = abs(fft(x)).^2;
Sy = abs(fft(y)).^2;

%% Thikonov regularization
% smoothing operator
Gamma = toeplitz([2 -1 zeros(1,N-2)]);
Gamma(1,:) = [];
Gamma(end,:) = [];

Tlambda = 10.^(-4 : 0.2 : 10);
Nlambda = length(Tlambda);
Tab_MSE = zeros(1,Nlambda)+NaN;
Tab_J = zeros(1,Nlambda)+NaN;

figure
i=0;
for lambda = Tlambda 
    i = i+1;
    x_hat = (eye(N) + lambda * (Gamma'*Gamma))\y;
    
    Tab_MSE(i) = norm(x-x_hat)^2;
    Tab_J(i) = norm(y-x_hat)^2+lambda * norm(Gamma*x_hat);
    
    % as a function of lambda
    subplot(311)
    plot(vect_t,x,'g')
    hold on
    plot(vect_t,x_hat,'k')
    plot(vect_t,y)
    axis tight
    xlabel('t')
    legend('x','y','xhat')
    title(['lambda = ' num2str(Tlambda(i))])    
    hold off
    subplot(312)
    loglog(Tlambda,Tab_MSE);
    xlabel('$\lambda$','interpreter','latex')
    ylabel('$||x-\hat{x}_{\lambda}||_2^2$','interpreter','latex')
    xlim([Tlambda(1) Tlambda(end)]);
    subplot(313)
    loglog(Tlambda,Tab_J);
    xlabel('$\lambda$','interpreter','latex')
    ylabel('$J(\hat{x}_{\lambda})$','interpreter','latex')
    xlim([Tlambda(1) Tlambda(end)]);
    drawnow
    
end
close

% looking for optimal lambda
[minerr, i_opt] = min(Tab_MSE);
lambda_opt = Tlambda(i_opt);
x_hat_opt = (eye(N) + lambda_opt * (Gamma'*Gamma))\y;

% spectral analysis - periodogram
Sx_hat_opt = abs(fft(x_hat_opt)).^2;

%% plots
% temporal
figure
plot(vect_t,x,'g')
hold on
plot(vect_t,y)
plot(vect_t,x_hat_opt,'k')
hold off
legend('x','y','xhat')
axis tight

% spectral analysis
nfft = (-N/2:N/2-1)*fe/N;
figure
semilogy(nfft,fftshift(Sx),'g')
hold on
semilogy(nfft,fftshift(Sy))
semilogy(nfft,fftshift(Sx_hat_opt),'k')
% axis([nfft(1) nfft(end) Sx(1) max(Sx)])
xlabel(['Frequency (f_e=' int2str(fe) ')'])
ylabel('Power spectral density')
legend('x','y','xhat')
axis tight
hold off

if isprint
    print([graphicpath 'fig_part4_denoising_PSD.png'],'-dpng','-r256')
end

% as a function of lambda
figure
subplot(311)
plot(vect_t,x,'g')
hold on
plot(vect_t,x_hat_opt,'k')
plot(vect_t,y)
axis tight
xlabel('t')
legend('x','y','xhat')
hold off
subplot(312)
loglog(Tlambda,Tab_MSE);
hold on
plot(lambda_opt,Tab_MSE(i_opt),'ro')
text(lambda_opt/10,minerr*10,['\lambda_{opt} = ' num2str(lambda_opt) ],'position',([0.61 0.8]),'color','r')
hold off
xlabel('$\lambda$','interpreter','latex')
ylabel('$||x-\hat{x}_{\lambda}||_2^2$','interpreter','latex')
subplot(313)
loglog(Tlambda,Tab_J);
hold on
plot(lambda_opt,Tab_J(i_opt),'ro')
hold off
xlabel('$\lambda$','interpreter','latex')
ylabel('$J(\hat{x}_{\lambda})$','interpreter','latex')

