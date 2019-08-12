function [y] = test_neumann(x)
%TEST_NEUMANN Summary of this function goes here
%   Detailed explanation goes here
tmp1=x(:,2)==0;
tmp2=x(:,1)==1;
tmp3=x(:,2)==1;
tmp4=x(:,1)==0;

y=zeros(length(x),1);
y(tmp1)=0;%x(tmp1,1);
y(tmp2)=1;%x(tmp2,2);
y(tmp3)=1;%x(tmp3,1);
y(tmp4)=0;%x(tmp4,2);
y=2*y;
end

