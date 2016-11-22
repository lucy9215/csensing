clear all
% M=8;
% N=12;
% Phi=PartHadamardMtx(M,N);
% 
% 
% 
% A=rgb2gray(imread('lena.jpg'));
% figure;
% imshow(A);
% 
% F=fft2(A);
% Fim=log10(abs(F));
% figure;
% imshow(Fim,[]);
% 
% S=fftshift(abs(F));
% S_norm=log10(S);
% figure;
% imshow(S_norm,[]);

n=16;
N=n*n;
x=[-ones(n-1,n); zeros(1,n)];
xl=reshape([-ones(n-1,n); zeros(1,n)],N,1);
Dv=[reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
  reshape([zeros(1,n); ones(n-1,n)],N,1)];

Dvv=spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
  reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);

Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
  reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);




