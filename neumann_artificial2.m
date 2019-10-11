function [y] = neumann_artificial2(x)
%TEST_NEUMANN Summary of this function goes here
%   Detailed explanation goes here

tmp2=x(:,1)==1 | x(:,2)==0 | x(:,2)==1 ;
y=zeros(length(x),1);
y(tmp2)=1;%cos(2*pi*x(tmp2,2))*10;
figure(1113);hold on;plot(y(tmp2),'k')
end

