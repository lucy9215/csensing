%压缩感知重构算法测试CS_Reconstuction_KtoPercentage.m
%   绘制参考文献中的Fig.2
%   参考文献：Joel A. Tropp and Anna C. Gilbert 
%   Signal Recovery From Random Measurements Via Orthogonal Matching
%   Pursuit，IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,
%   DECEMBER 2007.
%   Elapsed time is 1448.966882 seconds.(@20150418night)
clear all;close all;clc;
%% 参数配置初始化
CNT = 1000;%对于每组(K,M,N)，重复迭代次数
N = 256;%信号x的长度
Psi = eye(N);%x本身是稀疏的，定义稀疏矩阵为单位阵x=Psi*theta
M_set = [52,100,148,196,244];%测量值集合
Percentage = zeros(length(M_set),N);%存储恢复成功概率
%% 主循环，遍历每组(K,M,N)
tic
for mm = 1:length(M_set)
    M = M_set(mm);%本次测量值个数
    K_set = 1:5:ceil(M/2);%信号x的稀疏度K没必要全部遍历，每隔5测试一个就可以了
    PercentageM = zeros(1,length(K_set));%存储此测量值M下不同K的恢复成功概率
    for kk = 1:length(K_set)
       K = K_set(kk);%本次信号x的稀疏度K
       P = 0;
       for cnt = 1:CNT %每个观测值个数均运行CNT次
            Index_K = randperm(N);
            x = zeros(N,1);
            x(Index_K(1:K)) = 5*randn(K,1);%x为K稀疏的，且位置是随机的                
            Phi = randn(M,N);%测量矩阵为高斯矩阵
            A = Phi * Psi;%传感矩阵
            y = Phi * x;%得到观测向量y
            theta = CS_OMP(y,A,K);%恢复重构信号theta
            x_r = Psi * theta;% x=Psi * theta
            if norm(x_r-x)<1e-6%如果残差小于1e-6则认为恢复成功 %残差是向量差的2范数
                P = P + 1;
            end
       end
       PercentageM(kk) = P/CNT*100;%计算恢复概率
    end
    Percentage(mm,1:length(K_set)) = PercentageM;
end
toc
save KtoPercentage1000test %运行一次不容易，把变量全部存储下来
%% 绘图
S = ['-ks';'-ko';'-kd';'-kv';'-k*'];
figure;
for mm = 1:length(M_set)
    M = M_set(mm);
    K_set = 1:5:ceil(M/2);
    L_Kset = length(K_set);
    plot(K_set,Percentage(mm,1:L_Kset),S(mm,:));%绘出x的恢复信号
    hold on;
end
hold off;
xlim([0 125]);
legend('M=52','M=100','M=148','M=196','M=244');
xlabel('Sparsity level(K)');
ylabel('Percentage recovered');
title('Percentage of input signals recovered correctly(N=256)(Gaussian)');