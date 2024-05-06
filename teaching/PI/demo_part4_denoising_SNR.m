clear all
close all
clc

%% Generation of the unknown signal

% 1D input (unknown) signal
p = 7;
N = 2^p;
vect_t = (1:N)';
fe = 500;
f0 = 1;
x = cos(2*pi*f0/fe*vect_t+pi/2);


%% Thikonov regularization
% smoothing operator
Gamma = toeplitz([2 -1 zeros(1,N-2)]);
Gamma(1,:) = [];
Gamma(end,:) = [];


%% analysis as a function of SNR
N_MC = 20;
TSNRdB = [-15:2:50];
Tlambda = 10.^(-4 : 0.2 : 10);
Tlambda_opt = zeros(N_MC,length(TSNRdB));
Px     = 1/N*norm(x)^2;

Nlambda = length(Tlambda);
Tab_MSE = zeros(1,Nlambda);

for nsimu = 1:N_MC
    nsimu
    j = 0;
    for SNRdB = TSNRdB
        j=j+1;

        % noisy data
        sigma2 = 10^(-SNRdB/10)*Px;
        b = sqrt(sigma2)*randn(N,1);
        y = x+b;

        % reoslution
        i=0;
        for lambda = Tlambda 
            i = i+1;
            x_hat = (eye(N) + lambda * (Gamma'*Gamma))\y;
            Tab_MSE(i) = norm(x-x_hat)^2;
        end

        % looking for optimal lambda
        [minerr, i_opt] = min(Tab_MSE);
        Tlambda_opt(nsimu,j) = Tlambda(i_opt);

    end
end

figure;
semilogy(TSNRdB,mean(Tlambda_opt));
xlabel('SNR (dB)')
ylabel('$\lambda_{\mathrm{opt}}$','interpreter','latex')
graphicpath = '../slides/figures/';
print([graphicpath 'fig_part4_denoising_SNR.png'],'-dpng','-r256')

