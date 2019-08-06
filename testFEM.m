
[node,elem,bdFlag]=rect_mesh(1,1,100,100); % triangulace
% bdFlag(bdFlag>2)=1; % na strane 2 Neumann, na ostatnich Dirichlet

%% VSTUPNI PARAMETRY f, g_D, g_N, k
f=@(x)(-1+0*x(:,1)); % zatizeni (konstanta +1)
%f=@(x)(-5*(x(:,1)-0.5).^2-5*(x(:,2)-0.5).^2); % zatizeni (parabola)

g_D=@(x)(0*x(:,1)); % Dirichlet (konstantni na 4, linearni na 1 a 3)
%g_D=@(x)(0*x(:,1)+1); % Dirichlet (konstanta 1)

g_N=@(x)(0*x(:,1)); % Neumann (konstanta -1)

k=ones(length(elem),1); % materialove konstanty
%k(1:end/2)=10;

%% METODA KONECNYCH PRVKU
tStart=tic;
u=FEM(node,elem,bdFlag,k,f,g_D,g_N); % MKP
toc(tStart)

%% VYKRESLENI
figure;
trisurf(elem,node(:,1),node(:,2),u);
axis equal