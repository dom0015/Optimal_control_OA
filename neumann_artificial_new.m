function [y] = neumann_artificial_new(x)
%TEST_NEUMANN Summary of this function goes here
%   Detailed explanation goes here

tmp2=x(:,1)==1;
y=zeros(length(x),1);
y(tmp2)=sin(4*pi*x(tmp2,2))*10;
%figure(1113);hold on;plot(y(tmp2),'k')
end

