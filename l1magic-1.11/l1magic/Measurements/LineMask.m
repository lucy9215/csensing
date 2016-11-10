% LineMask.m
%
% Returns the indicator of the domain in 2D fourier space for the 
% specified line geometry.
% Usage :  [M,Mh,mi,mhi] = LineMask(L,N)
%
% Written by : Justin Romberg
% Created : 1/26/2004
% Revised : 12/2/2004

function [M,Mh,mi,mhi] = LineMask(L,N)


thc = linspace(0, pi-pi/L, L);
%thc = linspace(pi/(2*L), pi-pi/(2*L), L);

M = zeros(N);

% full mask
for ll = 1:L

	if ((thc(ll) <= pi/4) | (thc(ll) > 3*pi/4))
		yr = round(tan(thc(ll))*(-N/2+1:N/2-1))+N/2+1;
    	for nn = 1:N-1
      	M(yr(nn),nn+1) = 1;
      end
  else 
		xc = round(cot(thc(ll))*(-N/2+1:N/2-1))+N/2+1;
		for nn = 1:N-1
			M(nn+1,xc(nn)) = 1;
		end
	end

end


% upper half plane mask (not including origin)
Mh = zeros(N);
Mh = M;
Mh(N/2+2:N,:) = 0; % 取上半部分，下半部分赋值0
Mh(N/2+1,N/2+1:N) = 0; % 去掉了一根线（应该是避免重复）


M = ifftshift(M);%ifftshift 对换数据左右两边（我看结果是对换了上下和左右） 把星形变成了四个角一样的
mi = find(M);    %find(M) 找出M中非零元素的位置
Mh = ifftshift(Mh); % 也是对换了上下和左右，于是留下了下半部分非零
mhi = find(Mh);

