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

%% subsampling
% crude approach
qx = 2;
qy = 2;

Y_2D = X_2D(1:qx:end,1:qy:end);

mask_2D = zeros(size(X_2D));
mask_2D(1:qx:end,1:qy:end) = 1;

mask = mask_2D(:);
ind0 = find(mask==0);
H = sparse(eye(Nrow*Ncol));
H(ind0,:) = [];


% matrix formulation
X = X_2D(:);
Y = H*X;
Y_2D_bis = reshape(Y,Nrow/qx,Ncol/qy);

figure
subplot(1,2,1); imshow(Y_2D,[]); title('Downsampled w/ crude appraoch)');
colormap('gray')
axis square; axis tight
subplot(1,2,2); imshow(Y_2D_bis,[]); title('Downsampled w/ matrix product');
colormap('gray')



%% naive inversion

X_hat = H\Y;
X_hat_2D = reshape(X_hat,Nrow,Ncol);

figure
subplot(2,2,1); imshow(X_2D,[]); title('Original');
colormap('gray')
subplot(2,2,2); imshow(Y_2D,[]); title('Downsampled');
colormap('gray')
subplot(2,2,4); imshow(X_hat_2D,[]); title('Super-resolved');
colormap('gray')


%% naive inversion after noisy observation

SNRdB = 80;
Px = 1/(Nrow*Ncol)*norm(X)^2;
sigma2 = 10^(-SNRdB/10)*Px;
Y_obs = Y+randn(Nrow*Ncol/(qx*qy),1)*sqrt(sigma2);

X_hat = H\Y_obs;
X_hat_2D = reshape(X_hat,Nrow,Ncol);

figure
subplot(2,2,1); imshow(X_2D,[]); title('Original');
colormap('gray')
subplot(2,2,2); imshow(Y_2D,[]); title('Downsampled');
colormap('gray')
subplot(2,2,3); imshow(reshape(Y_obs,Nrow/qx,Ncol/qy),[]); title(['Downsampled + Noisy [' num2str(SNRdB) 'db]']);
colormap('gray')
subplot(2,2,4); imshow(X_hat_2D,[]); title('Super-resolved');
colormap('gray')


