clear all
close all
clc

% fichier d'origine (image)
X_2D_full = double(imread('image/barbara.png'));
step = 2;
X_2D = X_2D_full(1:step:200,(end-200+1):step:end);

[Nrow Ncol] = size(X_2D);

figure(1)
imagesc(X_2D)
colormap('gray')
set(gca,'xtick',[],'ytick',[])

pause

% subsampling
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
subplot(1,2,1); imagesc(Y_2D); title('Subsampled (crude)');
colormap('gray')
axis square; axis tight
subplot(1,2,2); imagesc(Y_2D_bis); title('Subsampled (matrix product)');
colormap('gray')
axis square; axis tight

% inversion after noisy observation
SNRdB = -10;
Px = 1/(Nrow*Ncol)*norm(X)^2;
sigma2 = 10^(-SNRdB/10)*Px;
Y_obs = Y;+randn(Nrow/2*Ncol/2,1)*sqrt(sigma2);

X_hat = H\Y_obs;
X_hat_2D = reshape(X_hat,Nrow,Ncol);

pause

figure
subplot(2,2,1); imagesc(X_2D); title('Original');
colormap('gray')
axis square; axis tight
subplot(2,2,2); imagesc(Y_2D); title('Decimated');
colormap('gray')
axis square; axis tight
subplot(2,2,3); imagesc(reshape(Y_obs,Nrow/2,Ncol/2)); title('Decimated+Noisy');
colormap('gray')
axis square; axis tight
subplot(2,2,4); imagesc(X_hat_2D); title('Reconstructed');
colormap('gray')
axis square; axis tight

pause

% regularization
k=0;
dl = 4;
d2 = -1;
m = Nrow
for i=1:Nrow
    for j=1:Nrow
        A = findneigh(i,j,Nrow);
        nnb = size(A,1) ;
        for h=1:nnb
            a(k+h)= ij2k( i , j , m ) ;
            b(k+h)= ij2k(A(h,1),A(h,2),m) ;
            if h==1
                c(k+h) = dl ;
            else
                c(k+h) = d2;
            end
        end
        k = k+nnb;
    end
end

[tmp ind] = find(b>m^2);
b(ind) = [];
a(ind) = [];
c(ind) = [];
G = sparse(a,b,c,m^2,m^2);
lambda = 0.1;
X_hat_r = (H'*H + lambda*G'*G)\(H'*Y_obs);

X_hat_2D_r = reshape(X_hat_r,Nrow,Ncol);
figure
subplot(2,2,1); imagesc(X_2D); title('Original');
colormap('gray')
axis square; axis tight
subplot(2,2,2); imagesc(Y_2D); title('Convolved');
colormap('gray')
axis square; axis tight
subplot(2,2,3); imagesc(reshape(Y_obs,Nrow,Ncol)); title('Convolved+Noisy');
colormap('gray')
axis square; axis tight
subplot(2,2,4); imagesc(X_hat_2D_r); title('Deconvolved');
colormap('gray')
axis square; axis tight
