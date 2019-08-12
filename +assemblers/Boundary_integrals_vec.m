function [f] = Boundary_integrals_vec(tri_grid,Boundary_edges,f)
%BOUNDARY_INTEGRALS Summary of this function goes here
%   Detailed explanation goes here
node=tri_grid.node;
edge=tri_grid.edge;
edge_int=edge(Boundary_edges,:);

node_midpoints=(node(edge_int(:,1),:)+node(edge_int(:,2),:))/2;

if ~isnumeric(f)
   f=f(node_midpoints); 
end

N=size(node,1);     % number of nodes
%% VECTORIZATION lengths
ve=node(edge_int(:,1),:)-node(edge_int(:,2),:);
area=sqrt(sum(ve.^2,2));

%% 2x2 local mass matrix [1/3 1/6;1/6 1/3]
ii=edge_int(:,1);
jj=edge_int(:,2);

sA=1/2*area.*f;
idx_i=[ii;jj];
vals=[sA;sA];

f=sparse(idx_i,ones(length(idx_i),1),vals,N,1);

end

