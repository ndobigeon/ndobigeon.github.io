clear all
close all
clc

% eigenvectors
U = orth(rand(4,4));
V = orth(rand(4,4));

% noise-free observation
% y belongs to the orthogonal of u_4
y = null(U(:,4)')*rand(3,1)*100;

% noisy observation
dy  = rand(4,1)*1e-1;
y_n = y+dy;

%% Invertible well-conditioned system

% forward matrix
Sigma = diag([1 1 1 1]);
H = U*Sigma*V';

% checking the singular values
sort(eig(H*H'),'descend')

% checking the condition number
cond(H)

% condition number
kappa = Sigma(1)/Sigma(end)

% LS inversion
x_LS = H\y

% LS inversion
x_LSn = H\y_n
dx = x_LSn-x_LS;

% sensitivity
norm(dx)/norm(x_LS)

pause

%% Invertible ill-conditioned system

% forward matrix
Sigma = diag([1 1 1 1e-5]);
H = U*Sigma*V';

% checking the singular values
sort(eig(H*H'),'descend');

% condition number
kappa = Sigma(1)/Sigma(end)

% checking the condition number
cond(H)

% LS inversion
x_LS = H\y

% LS inversion
x_LSn = H\y_n
dx = x_LS-x_LSn;

% sensitivity
norm(dx)/norm(x_LS)