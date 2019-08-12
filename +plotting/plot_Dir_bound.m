function [f] = plot_Dir_bound(tri_grid,u,u_d,R_r)
%PLOT_DIR_BOUND Summary of this function goes here
%   Detailed explanation goes here

bound_node=R_r*tri_grid.node;

if bound_node(1,1)==bound_node(end,1)
    idx_axis=2;
else
    idx_axis=1;
end



f=figure;
subplot(1,2,1);
hold on
plot(bound_node(:,idx_axis),R_r*u) % dirichlet condition vs result on Dir. boundary
plot(bound_node(:,idx_axis),u_d)
legend({'u','u_d'})
grid on
subplot(1,2,2);
plot(bound_node(:,idx_axis),abs(R_r*u-u_d))
set(gca,'yScale','log')
legend({'absolute_difference'})
grid on
end

