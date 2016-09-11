%% 验证OMP残差求解过程与Schmidt正交化的关系
M = 4;N = 10;
Phi = randn(M,N);
for nn = 1:N
    Phi(:,nn) = Phi(:,nn)/norm(Phi(:,nn));
end
b = randn(M,1);
e0 = b;%初始化残差为待稀疏信号b
%选第1列
c1 = Phi'*e0;%求出矩阵Phi每列与b的内积
[val1,pos1]=max(abs(c1));%找到内积中最大的列及其内积值
phit = [Phi(:,pos1)];%由所有选出的列组合的矩阵
Pphi = phit*(phit'*phit)^(-1)*phit';%正交投影变换矩阵
e1ompe = e0 - Pphi*e0;%OMP用上一次残差减去残差在phit列空间的正交投影
e1ompb = b - Pphi*b;%OMP用待稀疏信号b减去b在phit列空间的正交投影
x = Phi(:,pos1);%Schimidt正交化第一个向量
Px = x*(x'*x)^(-1)*x';
%实际上是b - Px*b
e1ompsmt = e0 - Px*b;
e1 = e1ompe;
norm(e1ompe-e1ompb)+norm(e1ompsmt-e1ompb)
%选第2列
c2 = Phi'*e1;%求出矩阵Phi每列与e1的内积
[val2,pos2]=max(abs(c2));%找到内积中最大的列及其内积值
phit = [Phi(:,pos1) Phi(:,pos2)];%由所有选出的列组合的矩阵
Pphi = phit*(phit'*phit)^(-1)*phit';%正交投影变换矩阵
e2ompe = e1 - Pphi*e1;%OMP用上一次残差减去残差在phit列空间的正交投影
e2ompb = b - Pphi*b;%OMP用待稀疏信号b减去b在phit列空间的正交投影
y = Phi(:,pos2) - Px*Phi(:,pos2);%Schimidt正交化第二个向量
Py = y*(y'*y)^(-1)*y';
%实际上是b - Px*b - Py*b
e2ompsmt = e1 - Py*b;%上一次残差减去b在第2列正交化所得z上的投影
e2 = e2ompe;
norm(e2ompe-e2ompb)+norm(e2ompsmt-e2ompb)
%选第3列
c3 = Phi'*e2;%求出矩阵Phi每列与e2的内积
[val3,pos3]=max(abs(c3));%找到内积中最大的列及其内积值
phit = [Phi(:,pos1) Phi(:,pos2) Phi(:,pos3)];%由所有选出的列组合的矩阵
Pphi = phit*(phit'*phit)^(-1)*phit';%正交投影变换矩阵
e3ompe = e2 - Pphi*e2;%OMP用上一次残差减去残差在phit列空间的正交投影
e3ompb = b - Pphi*b;%OMP用待稀疏信号b减去b在phit列空间的正交投影
z = Phi(:,pos3) - Px*Phi(:,pos3) - Py*Phi(:,pos3);%Schimidt正交化第三个向量
Pz = z*(z'*z)^(-1)*z';
%实际上是b - Px*b - Py*b - Pz*b
e3ompsmt = e2 - Pz*b;%上一次残差减去b在第3列正交化所得z上的投影
e3 = e3ompe;
norm(e3ompe-e3ompb)+norm(e3ompsmt-e3ompb)