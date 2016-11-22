% tveq_logbarrier.m
%
% 对应tv P_1 等式约束的logbarrier求解方法
% 
% Solve equality constrained TV minimization
% min TV(x)  s.t.  Ax=b.
%
% Recast as the SOCP
% min sum(t) s.t.  ||D_{ij}x||_2 <= t,  i,j=1,...,n
%                  Ax=b
% and use a log barrier algorithm.
%
% Usage:  xp = tveq_logbarrier(x0, A, At, b, lbtol, mu, slqtol, slqmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - Either a handle to a function that takes a N vector and returns a K 
%     vector , or a KxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes a K vector and returns an N vector.
%      If A is a KxN matrix, At is ignored.
%
% b - Kx1 vector of observations.
%
% lbtol - The log barrier algorithm terminates when the duality gap <= lbtol.
%         Also, the number of log barrier iterations is completely
%         determined by lbtol.
%         Default = 1e-3.   % 可忍受的误差，在此误差下停止迭代
%
% mu - Factor by which to increase the barrier constant at each iteration.
%      Default = 10.   迭代的步长
%
% slqtol - Tolerance for SYMMLQ; ignored if A is a matrix.
%     Default = 1e-8.   可忍受的误差 如果是矩阵则忽略此参数，因为矩阵不是largescal
%
% slqmaxiter - Maximum number of iterations for SYMMLQ; ignored
%     if A is a matrix.
%     Default = 200.    最大迭代次数
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function xp = tveq_logbarrier(x0, A, At, b, lbtol, mu, slqtol, slqmaxiter)  

% 判断是不是large scale模式，如果是largescale为1.这里如果A用的是fuction handle，那么判定为启用large scale模式
largescale = isa(A,'function_handle');  

% 默认参数的判断
if (nargin < 5), lbtol = 1e-3; end      % logbarrier停止迭代的误差
if (nargin < 6), mu = 10; end           % 步长默认10
if (nargin < 7), slqtol = 1e-8; end     % 可以忍受的误差
if (nargin < 8), slqmaxiter = 200; end  % 最大迭代次数

newtontol = lbtol;    % 停止迭代的误差
newtonmaxiter = 50;   % 最大迭代次数

N = length(x0);       % 传入的是一维向量，这里N=nxn
n = round(sqrt(N));   % 通过N计算出n

% create (sparse) differencing matrices for TV
% spdiags 取出非零元素
% [-ones(n-1,n); zeros(1,n)]是n*n的矩阵，n-1*n为-1，最后一行为0
% reshape以后按列取，变成n-1个-1，然后一个0重复...的一维向量
% [reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
%  reshape([zeros(1,n); ones(n-1,n)],N,1)]就成了两行向量
% 经过spdiags函数以后，Dv就成了对角线是-1从(1,1)-(63,63)，对角线左上的斜线是1的稀疏矩阵(1,2)(63,64)
% Dh是对角线为-1 距离为16的地方(1,17)开始是1的稀疏矩阵

Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
  reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
  reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);

% starting point --- make sure that it is feasible
if (largescale)
  if (norm(A(x0)-b)/norm(b) > slqtol) % 对传入的x0投影，减去观测值得到残差，求残差和观测值的2范数比例 大于可忍受的
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    AAt = @(z) A(At(z));
    [w,cgres] = cgsolve(AAt, b, slqtol, slqmaxiter, 0); % 这里用到了函数cgsolve
    if (cgres > 1/2)                   % 检查是否满足求解条件
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = At(w); % 求解后的w反投影为x0
  end
else    %　小规模下是否满足求解条件
  if (norm(A*x0-b)/norm(b) > slqtol)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    opts.POSDEF = true; opts.SYM = true;
    [w, hcond] = linsolve(A*A', b, opts);
    if (hcond < 1e-14)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = A'*w; % 对解出的w反投影为x0
  end  
end
x = x0;   % 为x赋初始值
Dhx = Dh*x;  Dvx = Dv*x; %　求gradient　horizontal and vertical
t = (0.95)*sqrt(Dhx.^2 + Dvx.^2) + (0.1)*max(sqrt(Dhx.^2 + Dvx.^2)); % 加权方法算出来的gradient值

% choose initial value of tau so that the duality(二元) gap after the first
% step will be about the origial TV
tau = N/sum(sqrt(Dhx.^2+Dvx.^2));   % 这个tau不知道是什么

lbiter = ceil((log(N)-log(lbtol)-log(tau))/log(mu));   %计算logbarrier的迭代次数
disp(sprintf('Number of log barrier iterations = %d\n', lbiter));
totaliter = 0;
for ii = 1:lbiter  % 进行logbarrier迭代
  
  [xp, tp, ntiter] = tveq_newton(x, t, A, At, b, tau, newtontol, newtonmaxiter, slqtol, slqmaxiter); % 进行迭代 返回值中有迭代次数
  totaliter = totaliter + ntiter;  %总迭代次数由牛顿迭代的次数累加
  
  tvxp = sum(sqrt((Dh*xp).^2 + (Dv*xp).^2)); % 总的gradient值 
  disp(sprintf('\nLog barrier iter = %d, TV = %.3f, functional = %8.3f, tau = %8.3e, total newton iter = %d\n', ...
    ii, tvxp, sum(tp), tau, totaliter));  % 这里的functional不知道是什么
  
  x = xp;
  t = tp;
  
  tau = mu*tau; %　新的tau是步长mu*tau
  
end
                   