function [R] = Restriction_other(tri_grid,node_idx)
%RESTRICTION Summary of this function goes here
%   Detailed explanation goes here
n=length(tri_grid.node);

if length(node_idx)~=n
    tmp=false(n,1);
    tmp(node_idx)=true;
    f_idx=~tmp;
else
    f_idx=~node_idx;
end
tmp=speye(n,n);
R=tmp(:,f_idx);
end

