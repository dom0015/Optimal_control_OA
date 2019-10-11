function [res]=get_neumann_triangles(tri_grid,R_m,diff_x)


n=length(tri_grid.node);
m=length(tri_grid.elem);

B=sparse([1:m 1:m 1:m],[tri_grid.elem(:,1) tri_grid.elem(:,2) tri_grid.elem(:,3)],ones(3*m,1),m,n);

R_m_tri=sum(B*R_m',2)==2;
res=diff_x(R_m_tri);
end