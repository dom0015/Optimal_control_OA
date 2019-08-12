function [f] = plot_res(u,tri_grid)
%PLOT_RES Summary of this function goes here
%   Detailed explanation goes here
f=figure;
h=trisurf(tri_grid.elem,tri_grid.node(:,1),tri_grid.node(:,2),full(u));

h.EdgeColor = 'none';
colormap jet(1000)
axis equal
view(0,90)
colorbar
end

