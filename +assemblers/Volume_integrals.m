function [A,M,f] = Volume_integrals( tri_grid, mat_stifness, mat_mass,f)
%FEM solution of -div(k*grad(p))=f
%   node ... coordinates of te nodes
%   elem ... vertices of the elements
%   dbFlag ... types of the sides
%               0 ... non-boundary side
%               1 ... Dirichlet boundary condition
%               2 ... Neumann boundary condition
node=tri_grid.node;
elem=tri_grid.elem;
k_stiff=mat_stifness;
k_mass=mat_mass;

N=size(node,1);     % number of nodes
NT=size(elem,1);    % number of elements

if length(k_stiff)==1
    k_stiff=ones(NT,1)*k_stiff;
end

if length(k_mass)==1
    k_mass=ones(NT,1)*k_mass;
end

%% VECTORIZATION (EDGES + AREA)
ve=zeros(NT,2,3);
ve(:,:,1)=node(elem(:,3),:)-node(elem(:,2),:);
ve(:,:,2)=node(elem(:,1),:)-node(elem(:,3),:);
ve(:,:,3)=node(elem(:,2),:)-node(elem(:,1),:);

area=0.5*abs(-ve(:,1,3).*ve(:,2,2)+ve(:,2,3).*ve(:,1,2));

%% USING THE SYMMETRY + PER PARTES
ii=zeros(NT,1,'double'); 
jj=zeros(NT,1,'double');
A=sparse(N,N);
M=sparse(N,N);
for i=1:3
    ii(:)=elem(:,i);
    sA=k_stiff.*dot(ve(:,:,i),ve(:,:,i),2)./(4*area);
    A=A+sparse(ii,ii,sA,N,N);
    sM=k_mass.*(2*area)/12;
    M=M+sparse(ii,ii,sM,N,N);
end
for i=1:2
    for j=i+1:3
        ii(:)=elem(:,i);
        jj(:)=elem(:,j);
        sA=k_stiff.*dot(ve(:,:,i),ve(:,:,j),2)./(4*area);
        A=A+sparse(ii,jj,sA,N,N);
        A=A+sparse(jj,ii,sA,N,N);
        sM=k_mass.*(2*area)/24;
        M=M+sparse(ii,jj,sM,N,N);
        M=M+sparse(jj,ii,sM,N,N);
    end
end

%% volume force/RIGHT HAND SIDE - 3-POINTS QUADRATURE
%tic;
mid1=(node(elem(:,2),:)+node(elem(:,3),:))/2; % midpoints
mid2=(node(elem(:,1),:)+node(elem(:,3),:))/2; % midpoints
mid3=(node(elem(:,1),:)+node(elem(:,2),:))/2; % midpoints

bt1=area.*(f(mid2)+f(mid3))/6;
bt2=area.*(f(mid1)+f(mid3))/6;
bt3=area.*(f(mid1)+f(mid2))/6;

f=accumarray(elem(:),[bt1;bt2;bt3],[N 1]);

end

