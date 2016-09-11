t = 0.1:0.1:10;
f1= 10;
f2= 15;

zer=ones(1,25);
one=zeros(1,25);
mask=[one,zer,one,zer];
s = sin(2*pi*f1*t).*[one,zer,zer,zer]+sin(2*pi*f1*t).*[zer,one,zer,zer]+sin(2*pi*f1*t).*[zer,zer,one,zer];


[sa,sd]=dwt(s,'haar');
m=[sa;sd];
% x=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.2*cos(2*pi*f3*Ts*ts)+0.9*cos(2*pi*f4*Ts*ts)+noise(Ts);  %  完整信号,由4个信号叠加而来
