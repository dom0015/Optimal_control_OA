function [A] = Boundary_integrals_mat(tri_grid,Boundary_edges,mat_const)
%BOUNDARY_INTEGRALS Summary of this function goes here
%   Detailed explanation goes here
node=tri_grid.node;
edge=tri_grid.edge;
edge_int=edge(Boundary_edges,:);


N=size(node,1);     % number of nodes
%% VECTORIZATION lengths
ve=node(edge_int(:,1),:)-node(edge_int(:,2),:);
area=sqrt(sum(ve.^2,2));

%% 2x2 local mass matrix [1/3 1/6;1/6 1/3]
ii=edge_int(:,1);
jj=edge_int(:,2);

sA=1/3*mat_const.*area;
idx_i=[ii;ii;jj;jj];
idx_j=[ii;jj;ii;jj];
vals=[sA;sA/2;sA/2;sA;];

A=sparse(idx_i,idx_j,vals,N,N);
end

