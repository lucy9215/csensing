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
X = phantom(n); % phantom的图像 大小n*n
x = X(:);       % 把phantom展开成一维向量（n*n）*1

% number of radial lines in the Fourier domain
L = 22;

% Fourier samples we are given
[M,Mh,mi,mhi] = LineMask(L,n);
OMEGA = mhi;    % mhi其实只有一半
A = @(z) A_fhp(z, OMEGA);      % 函数已输入OMEGA    调用A(x) 相当于调用A_fhp(x,OMEGA)
At = @(z) At_fhp(z, OMEGA, n); % 函数已输入OMEGA和n 调用At(y) 相当于调用A_fhp(y,OMEGA,n)

% measurements
y = A(x); % 测量值 OMEGA个img和OMEGA个real以及一个均值，是mask以后的频域

% min l2 reconstruction (backprojection)
xbp = At(y); % minimal energy reconstruction initial guess 初始猜测的最小能量恢复 是用mask以后的频域恢复的图像向量
Xbp = reshape(xbp,n,n); % 最小能量恢复的图 reshape把向量组织成图



% recovery
tic
tvI = sum(sum(sqrt([diff(X,1,2) zeros(n,1)].^2 + [diff(X,1,1); zeros(1,n)].^2 ))); % 最开始的TV值
% TV的值是每个像素sqrt(横向的‘相邻元素的差值’^2+纵向‘相邻元素的差值’^2)，总nxn个像素的gradient总和
disp(sprintf('Original TV = %8.3f', tvI));
xp = tveq_logbarrier(xbp, A, At, y, 1e-1, 2, 1e-8, 600); % 进行优化
Xtv = reshape(xp, n, n); % 范数优化输出的图像
toc

figure;
subplot(1,2,1);
title('minimal energy reconstruction');
imshow(Xbp);
subplot(1,2,2);
title('l1-magic reconstruction');
imshow(Xtv);
