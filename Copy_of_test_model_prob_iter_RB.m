%% Discretization size and problem setting
nx=60; % discretization size x axis
ny=60; % y axis
Lx=1;   % size of the domain x axis
Ly=1;   % y axis
smoothing_boundary = false; % deforms grid to be smoother around the boundary
% this makes neumann boundary conditions much
% more precise

% setting for the problem with solution: u(x,y)= x*x+y*y (dirichlet match
% the solution), therefore for low beta it match with Prescribed DB

f=@(x)x(:,1)*0+1;            % rhs of PDE
Neumann_boundary=@(x)0*x(:,1);       % neuman boundary conditions (g(x))
Dirichlet_boundary={@(x)x(:,1)*0, @(x)x(:,1)*0, @(x)x(:,1)*0}; % Dirichlet boundary to match
% Given neumann is on top,left and bottom side
% Given Dirichlet to optimize is on the left side (0.2,0.8) part of side
b_Dir={4,[0 1]};
b_Neu_known={4,[0 1]};
b_Neu_unknown={3,[0 1]
    1,[0 1]
    2,[0 1]};

% other parameters
sigma=1;


%% Assemble all matrices and vector of the problem
[M_r,M_m,K,R_r,R_m,R_b,f_vec,g_vec,u_d,tri_grid] = assemblers.Assembly_all(nx,ny,Lx,Ly,...
    f,Neumann_boundary,Dirichlet_boundary,sigma,b_Dir,b_Neu_known,b_Neu_unknown,smoothing_boundary);


%% add artificially computed u_d
u_d = compute_artificial_u_d_2(nx,ny,smoothing_boundary);

%% Assembling 3x3 block matrix
n_u=length(tri_grid.node);
n_v=size(R_m,1);

M_r_pruh=R_r'*M_r*R_r;
N=R_m'*M_m;


%beta=1e-9;

neuman_norm=[];
dirichlet_norm=[];
min_eig=[];
eig_comp_real=[];
eig_comp_imag=[];
kkk=0;
iter_count=[];
res_res=[];

for beta=10.^(-(1:9))
    kkk=kkk+1;
    
    A_3x3= [M_r_pruh            K             sparse(n_u,n_v)
        K              sparse(n_u,n_u)   -N
        sparse(n_v,n_u)    -N'            beta*M_m];
    b_3x3=[R_r'*M_r*u_d
        f_vec+R_b'*g_vec
        sparse(n_v,1)];
    
    % solution
    res=A_3x3\b_3x3;
    
    % iterative
    %M_m_diag_inv=diag(1./diag(M_m));
    M_m_pruh=R_m'*M_m*R_m;
    
    A_2x2=[K         (-1/beta)*M_m_pruh
        M_r_pruh                  K];
    B_2x2=[K         (-1/beta)*M_m_pruh
        M_r_pruh  K+M_r_pruh+(1/beta)*M_m_pruh];
    b_2x2=[f_vec+R_b'*g_vec
        R_r'*M_r*u_d];
    ddim=size(b_2x2);
    
    % iterative solution fgmres/gmres)
    restart_iter=100; tol = 1e-6;  maxit = 100;
    [x,iter,resvec]=itersolvers.fgmres(A_2x2,b_2x2,restart_iter,tol,maxit,B_2x2);
    % [x,flag,relres,iter,resvec]=gmres(A_2x2,b_2x2,restart_iter,tol,maxit,B_2x2);
    % iter
    % resvec
    
    % spectrum estimation for discretization less than 20 it makes full
    % eigenvalue decomp. for greater, it computes smallest 20 eigenvalues
    
    if nx<=20
        PA=B_2x2\A_2x2;
        [U,D]=eig(full(PA));
        EI=diag(D);
    else
        EI=eigs(A_2x2,B_2x2,nx*2+3,'smallestabs');
    end
    
    % plotting eigenvalues/lowest 20
    figure
    [tmp,idx]=sort(real(EI));
    hold on
    plot(imag(EI(idx)),'r.','MarkerSize',15)
    plot(tmp,'k.-','MarkerSize',10)
    grid on
    %set(gca,'xScale','log')
    legend({'Imaginary part of Eig(PA)','Real part of Eig(PA)'})
    
    
    %%  extracting variables
    u=res(1:n_u);
    w=res((n_u+1):2*n_u);
    v=res((2*n_u+1):end);
    
    %% plotting results
    plotting.plot_res(u,tri_grid); % values of u
    [f,diff_x,diff_y]=plotting.plot_res_grad(u,tri_grid); % gradient of u
    plotting.plot_Dir_bound(tri_grid,u,u_d,R_r); % difference on Dir. Boundary
    [res]=get_neumann_triangles(tri_grid,R_m,diff_x);
    figure(1113);hold on; plot(res)
    neuman_norm(kkk)=norm(res-1);
    dirichlet_norm(kkk)=norm(R_r*u-u_d);
    min_eig(kkk)=tmp(1);
    [~,idx]=max(imag(EI));
    eig_comp_real(kkk)=real(EI(idx));
    eig_comp_imag(kkk)=imag(EI(idx));
    iter_count(kkk)=iter;
    res_res(kkk)=resvec(end);
end
