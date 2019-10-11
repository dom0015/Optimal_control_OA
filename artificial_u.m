function [u,u_d] = artificial_u (nx,ny,Lx,Ly,smothing_boundary,Neumann_boundary,f,sigma)

%% Discretization size and broblem setting

Dirichlet_boundary=@(x)x(:,1)*0; % Dirichlet boundary to match
% Given neumann is on top,left and bottom side
% Given Dirichlet to optimize is on the left side (0.2,0.8) part of side
b_Dir={4,[0 1]};
b_Neu_known={1,[0 1]
    2,[0 1]
    3,[0 1]
    4,[0 1]};
b_Neu_unknown={1,[0 1]};

% other parameters


%% Assemble all matrices and vector of the problem
[~,~,K,R_r,~,R_b,f_vec,g_vec,~,~] = assemblers.Assembly_all(nx,ny,Lx,Ly,...
    f,Neumann_boundary,Dirichlet_boundary,sigma,b_Dir,b_Neu_known,b_Neu_unknown,smothing_boundary);

u=K\(f_vec+R_b'*g_vec);

%% plotting results
%plotting.plot_res(u,tri_grid); % values of u
%plotting.plot_res_grad(u,tri_grid); % gradient of u

u_d=R_r*u;
end