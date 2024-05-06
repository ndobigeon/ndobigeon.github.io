clear all
close all
clc

N = 6;
M = 4;

U = orth(rand(N,N));
V = orth(rand(M,M));

Sigma = zeros(N,M);
D = min(M,N);

Tsigma = [ones(1,D)]; Tsigma(end) = 0;

for i=1:D
    Sigma(i,i) = Tsigma(i);
end

H = U*Sigma*V';
disp('rank H')
rank(H)

x = randn(M,1);

y = H*x;

xhat = H\y;

disp('normalized error')
norm(x-xhat)/norm(x)

disp('cond H')
cond(H)

disp('residual error')
norm(y-H*xhat)/norm(y)

Hdag = (H'*H)^(-1)*H';
Hdag*H
H*Hdag