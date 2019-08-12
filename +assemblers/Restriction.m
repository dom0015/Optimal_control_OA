function [R] = Restriction(tri_grid,node_idx)
%RESTRICTION Summary of this function goes here
%   Detailed explanation goes here
n=length(tri_grid.node);
tmp=speye(n,n);
R=tmp(:,node_idx);
end

