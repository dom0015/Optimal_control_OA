function u_d = compute_artificial_u_d_1(nx,ny,smothing_boundary)

%% Discretization size and broblem setting
Lx=1;   % size of the domain x axis
Ly=1;   % y axis

% setting for the problem with solution: u(x,y)= x*x+y*y (dirichlet match
% the solution), therefore for low beta it match with Prescribed DB

f=@(x)x(:,1)*0+1;            % rhs of PDE
Neumann_boundary=@(x)neumann_artificial1(x);       % neuman boundary conditions (g(x))
Dirichlet_boundary={@(x)x(:,1)*0, @(x)x(:,1)*0, @(x)x(:,1)*0}; % Dirichlet boundary to match
% Given neumann is on top,left and bottom side
% Given Dirichlet to optimize is on the left side (0.2,0.8) part of side
b_Dir={3,[0 1]
    4,[0 1]
    1,[0 1]};
b_Neu_known={3,[0 1]
    4,[0 1]
    1,[0 1]
    2,[0 1]};
b_Neu_unknown={2,[0 1]};

% other parameters
sigma=1;

%% Assemble all matrices and vector of the problem
[~,~,K,R_r,~,R_b,f_vec,g_vec,~,tri_grid] = assemblers.Assembly_all(nx,ny,Lx,Ly,...
    f,Neumann_boundary,Dirichlet_boundary,sigma,b_Dir,b_Neu_known,b_Neu_unknown,smothing_boundary);

u=K\(f_vec+R_b'*g_vec);

%% plotting results
%plotting.plot_res(u,tri_grid); % values of u
%plotting.plot_res_grad(u,tri_grid); % gradient of u

u_d=R_r*u;
end