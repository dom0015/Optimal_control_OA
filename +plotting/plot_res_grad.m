function [f,diff_x,diff_y] = plot_res_grad(u,tri_grid)
%PLOT_RES Summary of this function goes here
%   Detailed explanation goes here
u=full(u);
x=tri_grid.node(:,1);
x=x(tri_grid.elem);
y=tri_grid.node(:,2);
y=y(tri_grid.elem);
z=u(tri_grid.elem);

tmp=(x(:,1).*y(:,2)-x(:,2).*y(:,1)-x(:,1).*y(:,3)+x(:,3).*y(:,1)+x(:,2).*y(:,3)-x(:,3).*y(:,2));
diff_x=-(y(:,1).*z(:,2)-y(:,2).*z(:,1)-y(:,1).*z(:,3)+y(:,3).*z(:,1)+y(:,2).*z(:,3)-y(:,3).*z(:,2))./tmp;
diff_y=(x(:,1).*z(:,2)-x(:,2).*z(:,1)-x(:,1).*z(:,3)+x(:,3).*z(:,1)+x(:,2).*z(:,3)-x(:,3).*z(:,2))./tmp;


f=figure;
subplot(1,2,1)
p=tri_grid.node';
t=tri_grid.elem';
x=p(1,:);
y=p(2,:);
P=[x(t(:));y(t(:))];
T=reshape(1:size(P,2),[3 size(P,2)/3]);
% create random u for testing
u=diff_x;
tmp=[u';u';u'];
h=trisurf(T',P(1,:),P(2,:),tmp(:));
h.EdgeColor = 'none';
colormap jet(1000)
axis equal
view(0,90)
colorbar

subplot(1,2,2)
u=diff_y;
tmp=[u';u';u'];
h=trisurf(T',P(1,:),P(2,:),tmp(:));
h.EdgeColor = 'none';
colormap jet(1000)
axis equal
view(0,90)
colorbar
end

