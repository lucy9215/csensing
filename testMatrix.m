clear all
M=8;
N=12;
Phi=PartHadamardMtx(M,N);



A=rgb2gray(imread('lena.jpg'));
figure;
imshow(A);

F=fft2(A);
Fim=log10(abs(F));
figure;
imshow(Fim,[]);

S=fftshift(abs(F));
S_norm=log10(S);
figure;
imshow(S_norm,[]);