function y=RIPText(x,t)  
% calculete the kth Restricted Isometry Constant of measurement matrix x.  
% Modified on March 13th, 2012.  
% Reference: Wei Dai, Olgica Milenkovic. Subspace Pursuit for Compressive Sensing  
% Signal Reconstruction. IEEE Trans. ona Information Theory. vol.55, no.5,  
% 2230-2249. May, 2009  
%calculate the Restricted Isomentry Constant of matrix constructed by  
%random n columns of the original matrix x.  
%created by Li Bo in March 16th, 2011  
%Harbin Institute of Technolgy  
n=size(x,2);%the numbers of column of original matrix x  
Num=1:n;%create a row of numbers from one to n with iterval 1  
Com=combntns(Num,t);% all the combination of the selected t columns from total n columns  
ma=size(Com,1);% the number of combinations  
delta=zeros(1,ma);% a vector used to store the restricted isometry constant candidate  
for i=1:ma  
    b=x(:,Com(i,:));% new matrix constructed by the t selected columns   
    e=eig(b'*b-eye(t));%calculate all the eigen values of matrix b'*b-eye(t) # b的特征值numda,这里得出numda-1
    delta(i)=(max(abs(e)));%select the largest absolute eigen value as restricted isomentry candidate  
end  
y=max(delta);% select the largest one as restricted isometry constant  