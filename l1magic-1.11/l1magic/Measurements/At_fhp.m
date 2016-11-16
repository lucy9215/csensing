% At_fhp.m
%
% Adjoint of At_fhp (2D Fourier half plane measurements).
%
% Usage: x = At_fhp(b, OMEGA, n)
%
% b - K vector = [mean; real part(OMEGA); imag part(OMEGA)]
%
% OMEGA - K/2-1 vector denoting which Fourier coefficients to use
%         (the real and imag parts of each freq are kept).
%
% n - Image is nxn pixels
%
% x - N vector
%
% Written by: Justin Romberg, Caltech
% Created: October 2005
% Email: jrom@acm.caltech.edu
%

function x = At_fhp(y, OMEGA, n)

K = length(y);

fx = zeros(n,n);  % 初始化频域
fx(1,1) = y(1);   % 初始化均值
fx(OMEGA) = sqrt(2)*(y(2:(K+1)/2) + i*y((K+3)/2:K)); % 在频域的OMEGA处赋值 % 这里的y其实是一维向量，所以从中取出对应的OMEGA位置的向量
x = reshape(real(n*ifft2(fx)), n*n, 1);    % 对此
