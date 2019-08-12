function [y] = neumann_artificial(x)
%TEST_NEUMANN Summary of this function goes here
%   Detailed explanation goes here

tmp2=x(:,1)==1 & x(:,2)>=0.15 & x(:,2)<=0.85 ;
y=zeros(length(x),1);
y(tmp2)=1;%x(tmp2,2);

end

