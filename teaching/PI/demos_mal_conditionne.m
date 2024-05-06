clear all
close all
clc

H = [15  5  6  5 ;...
      7  12  7  5 ;...
      1  2 14  3 ;...
      1  2  4 15 ];
eig(H*H')
D = sort(eig(H*H'),'descend');
C = D(1)/D(end)

y = [31 31 20 22]';

x = H\y

% perturbation
dy = [1 -1 1 -1]'*1e-1;
vary = norm(dy)

% sensibilité
x2 = H\(y+dy)
dx = x2-x;
varx = norm(dx)

pause

clear all

% mal conditionné
H = [10  7  8  7 ;...
      7  5  6  5 ;...
      8  6 10  9 ;...
      7  5  9 10 ];
eig(H*H')
D = sort(eig(H*H'),'descend');
C = D(1)/D(end)

y = [32 23 33 31]';
x = H\y

% perturbation
dy = [1 -1 1 -1]'*1e-1;
vary = norm(dy)

% sensibilité
x2 = H\(y+dy)
dx = x2-x;
varx = norm(dx)