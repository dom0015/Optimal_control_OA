function [M_r,M_m,K,R_r,R_m,R_b,f_vec,g_vec,u_d,tri_grid] = Assembly_all...
    (nx,ny,Lx,Ly,f,Neumann_boundary,Dirichlet_boundary,sigma,b_Dir,b_Neu_known,b_Neu_unknown,smothing_boundary)
%ASSEMBLY_ALL Summary of this function goes here
%   Detailed explanation goes here
%% Discretization size and broblem setting
mat_stifness=1;
mat_mass=sigma;


%% triangular discretization on rectangle
% trigrid is structure of node,elem,edge
% bound_idx is structure od node and edge idx on boundary
% (1-bottom,2-right,3-top,4-left))
[ tri_grid, bound_idx ] = grid.rect_mesh( Ly, Lx, ny, nx );


%% changing grid------------------------------------
% making grid more dense around the borders
if smothing_boundary
    tri_grid.node(:,1)=betacdf(tri_grid.node(:,1)/Lx,4,4)*Lx;
    tri_grid.node(:,2)=betacdf(tri_grid.node(:,2)/Ly,4,4)*Ly;
end
%----------------------------------------------------

%% SELECT boundary for Dirichlet/Neumann conditions
%b_Dir={4,[0.2 0.8]}; % nx2 cell array frist column - which side, second which part of the side scaled to <0,1>;
b_Dir_val= Dirichlet_boundary; % cell array of vector funnctions working on nx2 (x,y) inputs
[D_bound_set] = grid.D_boundary_setting(b_Dir,b_Dir_val,tri_grid,bound_idx);
% D_bound_set is structure of Dir node mask and their values

% Neumann known conditions
% b_Neu_known={3,[0 1]
%     4,[0 1]
%     1,[0 1]};
%b_Neu_val_known={Neumann_boundary;Neumann_boundary;Neumann_boundary};
%[N_bound_set] = grid.N_boundary_setting(b_Neu_known,b_Neu_val_known,tri_grid,bound_idx);

% Neumann unknown boundary
%b_Neu_unknown={2,[0 1]};

% creating boundary integration selections
[Dir_bound] = grid.Select_boundary(b_Dir,tri_grid,bound_idx);
[Neu_known_bound] = grid.Select_boundary(b_Neu_known,tri_grid,bound_idx);
[Neu_unknown_bound] = grid.Select_boundary(b_Neu_unknown,tri_grid,bound_idx);

%% Assembling matrices and rhss
% all originaly assembled matrices/rhs, are for all dof (size of nodes) -
% at default prolongated

% M/rhs coming from area integrations matrices are full sizes (as nodes)
[K_A,K_M,f_vec] = assemblers.Volume_integrals( tri_grid, mat_stifness, mat_mass,f);

% M/rhs coming from boundary integrations matrices are full sizes (as nodes)
M_r=assemblers.Boundary_integrals_mat(tri_grid,Dir_bound,1);
M_m=assemblers.Boundary_integrals_mat(tri_grid,Neu_unknown_bound,1);
g_vec=assemblers.Boundary_integrals_vec(tri_grid,Neu_known_bound,Neumann_boundary);

% restrictions to some sets of nodes
Neu_unknown_node_idx=tri_grid.edge(Neu_unknown_bound,:);
Neu_unknown_node_idx=unique(Neu_unknown_node_idx(:));
Neu_known_node_idx=tri_grid.edge(Neu_known_bound,:);
Neu_known_node_idx=unique(Neu_known_node_idx(:));
[R_m] = assemblers.Restriction(tri_grid,Neu_unknown_node_idx);
R_m=R_m';
[R_b] = assemblers.Restriction(tri_grid,Neu_known_node_idx);
R_b=R_b';
[R_r] = assemblers.Restriction(tri_grid,D_bound_set.node_flg);
R_r=R_r';

% building basic blocks
K=K_A+K_M;
M_r=R_r*M_r*R_r';
M_m=R_m*M_m*R_m';
u_d=R_r*D_bound_set.node_val;
g_vec=R_b*g_vec;
f_vec=f_vec;
end

