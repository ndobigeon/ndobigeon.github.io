clear all
close all
clc

%% fichier d'origine (image)
X_2D_full = double(imread('image/barbara.png'));
step = 2;
X_2D = X_2D_full(1:step:200,(end-200+1):step:end);

[Nrow Ncol] = size(X_2D);

figure(1)
imshow(X_2D,[])
colormap('gray')
set(gca,'xtick',[],'ytick',[])

%% gaussian blurring
figure
sigma = 1.8;
n_ker = 6;
psf  = fspecial('gaussian',n_ker+1,sigma);
subplot(2,2,1); surf(psf); title('PSF');
axis square; axis tight
subplot(2,2,2); imagesc(psf); title('PSF');
axis square; axis tight
       
Y_2D = imfilter(X_2D,psf);

%figure
subplot(2,2,3); imshow(X_2D,[]); title('Original');
colormap('gray')
subplot(2,2,4); imshow(Y_2D,[]); title('Convolved - Gaussian blur');
colormap('gray')

%% motion blurring
figure
n_ker = 6;
theta = -45;
psf2  = fspecial('motion',n_ker+2,theta);
subplot(2,2,1); surf(psf2); title('PSF');
axis square; axis tight
subplot(2,2,2); imagesc(psf2); title('PSF');
axis square; axis tight
       
Y_2D_motion = imfilter(X_2D,psf2);

subplot(2,2,3); imshow(X_2D,[]); title('Original');
colormap('gray')
subplot(2,2,4); imshow(Y_2D_motion,[]); title('Convolved - Motion blur');
colormap('gray')

%% matrix formulation
H = convmtx2(psf,Nrow,Ncol);
X = X_2D(:);

[N, M] = size(H);
masq1 = zeros(Nrow,Ncol) + 1;
masq2 = zeros(size(psf)+[Nrow Ncol]-1);

masq2(n_ker/2+(1:Nrow),n_ker/2+(1:Ncol)) = masq1;
ind = find(reshape(masq2,1,length(masq2(:)))~=0);
H = H(ind,:);  
Y = H*X;
Y_2D_bis = reshape(Y,Nrow,Ncol);

figure
subplot(1,2,1); imshow(Y_2D,[]); title('Convolved w/ imfilter');
colormap('gray')
subplot(1,2,2); imshow(Y_2D_bis,[]); title('Convolved w/ matrix product');
colormap('gray')

%% naive inversion

X_hat = H\Y;
X_hat_2D = reshape(X_hat,Nrow,Ncol);

figure
subplot(2,2,1); imshow(X_2D,[]); title('Original');
colormap('gray')
subplot(2,2,2); imshow(Y_2D,[]); title('Convolved');
colormap('gray')
subplot(2,2,4); imshow(X_hat_2D,[]); title('Deconvolved');
colormap('gray')


%% naive inversion after noisy observation

SNRdB = 80;
Px = 1/(Nrow*Ncol)*norm(X)^2;
sigma2 = 10^(-SNRdB/10)*Px;
Y_obs = Y+randn(Nrow*Ncol,1)*sqrt(sigma2);

X_hat = H\Y_obs;
X_hat_2D = reshape(X_hat,Nrow,Ncol);

figure
subplot(2,2,1); imshow(X_2D,[]); title('Original');
colormap('gray')
subplot(2,2,2); imshow(Y_2D,[]); title('Convolved');
colormap('gray')
subplot(2,2,3); imshow(reshape(Y_obs,Nrow,Ncol),[]); title(['Convolved + Noisy [' num2str(SNRdB) 'db]']);
colormap('gray')
subplot(2,2,4); imshow(X_hat_2D,[]); title('Deconvolved');
colormap('gray')


