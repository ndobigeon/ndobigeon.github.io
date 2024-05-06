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

% pause
% 
% % gaussian blurring
% figure
% sigma = 0.8;
% n_ker = 6;
% psf  = fspecial('gaussian',n_ker+1,sigma);
% subplot(1,2,1); surf(psf); title('PSF');
% axis square; axis tight
% subplot(1,2,2); imagesc(psf); title('PSF');
% axis square; axis tight
%        
% Y_2D = imfilter(X_2D,psf);
% figure
% subplot(1,2,1); imagesc(X_2D); title('Original');
% colormap('gray')
% axis square; axis tight
% subplot(1,2,2); imagesc(Y_2D); title('Filtered (Gaussian)');
% colormap('gray')
% axis square; axis tight

% pause

% motion blurring
figure
% n_ker = 6;
% theta = -45;
% psf  = fspecial('motion',n_ker+2,theta);
sigma = 0.8;
n_ker = 6;
psf  = fspecial('gaussian',n_ker+1,sigma);
subplot(1,2,1); surf(psf); title('PSF');
axis square; axis tight
subplot(1,2,2); imagesc(psf); title('PSF');
axis square; axis tight
       
Y_2D_motion = imfilter(X_2D,psf);
figure
subplot(1,2,1); imagesc(X_2D); title('Original');
colormap('gray')
axis square; axis tight
subplot(1,2,2); imagesc(Y_2D_motion); title('Filtered (Motion)');
colormap('gray')
axis square; axis tight

pause

% matrix formulation
H = convmtx2(psf,Nrow,Ncol);
X = X_2D(:);

[N M] = size(H);
masq1 = zeros(Nrow,Ncol) + 1;
masq2 = zeros(size(psf)+[Nrow Ncol]-1);
% masq2(3:34,3:34) = masq1;
masq2(n_ker/2+(1:Nrow),n_ker/2+(1:Ncol)) = masq1;
ind = find(reshape(masq2,1,length(masq2(:)))~=0);
H = H(ind,:);  
Y = H*X;
Y_2D_bis = reshape(Y,Nrow,Ncol);



% inversion after noisy observation
SNRdB = 35;
Px = 1/(Nrow*Ncol)*norm(X)^2;
sigma2 = 10^(-SNRdB/10)*Px;
Y_obs = Y+randn(Nrow*Ncol,1)*sqrt(sigma2);

X_hat = H\Y_obs;
X_hat_2D = reshape(X_hat,Nrow,Ncol);


% figure
% subplot(2,2,1); imshow(X_2D,[]); title('True image','fontsize',14);
% colormap('gray')
% axis square; axis tight
% % subplot(2,2,2); imagesc(Y_2D); title('Convolved');
% % colormap('gray')
% % axis square; axis tight
% subplot(2,2,2); imshow(reshape(Y_obs,Nrow,Ncol),[]); title('Observed image','fontsize',14);;%imagesc(reshape(Y_obs,Nrow,Ncol)); title('Convolved+Noisy');
% colormap('gray')
% axis square; axis tight
% subplot(2,2,4); imshow(X_hat_2D,[]); title('Restored image','fontsize',14);
% colormap('gray')
% axis square; axis tight

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
lambda = 0.001;
X_hat_r = (H'*H + lambda*G'*G)\(H'*Y_obs);

X_hat_2D_r = reshape(X_hat_r,Nrow,Ncol);

figure
subplot(1,4,1); imshow(X_2D,[]); title('True image','fontsize',12);
colormap('gray')
axis square; axis tight
% subplot(2,2,2); imagesc(Y_2D); title('Convolved');
% colormap('gray')
% axis square; axis tight
subplot(1,4,2); imshow(reshape(Y_obs,Nrow,Ncol),[]); title('Observed image','fontsize',12);;%imagesc(reshape(Y_obs,Nrow,Ncol)); title('Convolved+Noisy');
colormap('gray')
axis square; axis tight
subplot(1,4,3); imshow(X_hat_2D,[]); title('ML estimate','fontsize',12);
colormap('gray')
axis square; axis tight
subplot(1,4,4); imshow(X_hat_2D_r,[]); title('MAP estimate','fontsize',12);
colormap('gray')
axis square; axis tight
print('fig_Bayes_restoration.png','-dpng','-r512')

