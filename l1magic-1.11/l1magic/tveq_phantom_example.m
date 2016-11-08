% tveq_phantom_example.m
%
% Phantom reconstruction from samples on 22 radial lines in the Fourier
% plane.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%


path(path, './Optimization');
path(path, './Measurements');
path(path, './Data');


% Phantom 
n = 256;
N = n*n;
X = phantom(n);
x = X(:);

% number of radial lines in the Fourier domain
L = 22;

% Fourier samples we are given
[M,Mh,mh,mhi] = LineMask(L,n);
OMEGA = mhi;
A = @(z) A_fhp(z, OMEGA);
At = @(z) At_fhp(z, OMEGA, n);

% measurements
y = A(x); % 测量值

% min l2 reconstruction (backprojection)
xbp = At(y); % minimal energy reconstruction initial guess 初始猜测的最小能量恢复
Xbp = reshape(xbp,n,n); % 最小能量恢复的图



% recovery
tic
tvI = sum(sum(sqrt([diff(X,1,2) zeros(n,1)].^2 + [diff(X,1,1); zeros(1,n)].^2 )));
disp(sprintf('Original TV = %8.3f', tvI));
xp = tveq_logbarrier(xbp, A, At, y, 1e-1, 2, 1e-8, 600);
Xtv = reshape(xp, n, n);
toc

figure;
subplot(1,2,1);
title('minimal energy reconstruction');
imshow(Xbp);
subplot(1,2,2);
title('l1-magic reconstruction');
imshow(Xtv);



