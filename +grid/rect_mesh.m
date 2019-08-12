function [ tri_grid, bound_idx ] = rect_mesh( height, width, h_elem, w_elem )
%% Creation mesh (nodes and elements) from rectangle positioned in (0,0)
% Input height and width of rectangle, number of elements in row - w_elem,
% number elements in column - h_elem.
% Output node matrix = position of verticles of mash
%        elem matrix = indexes of nodes creating elements
%        bdFlag matrix = for each element contain border-flag for each edge 

%% nodes matrix creation
x1 = linspace(0,width,w_elem+1);
n2 = h_elem+1;
x2 = linspace(0,height,h_elem+1);
[X1,X2] = meshgrid(x1,x2);
node = [X1(:) X2(:)];
%% elem matrix creation
m = 2*w_elem*h_elem;
elem = zeros(m,3);
idx = 1;
temp=reshape( repmat( 1:h_elem, 2,1 ), 1, [] );
for i=1:w_elem
    first=temp+(i-1)*n2;    % vector of first verticles of elements
    second=[temp(2:end) (h_elem+1)]+i*n2; % second verticles
    third=temp+1+(i-1)*n2;  % third verticles 
    third(1:2:end)=third(1:2:end)+n2;
    % together make elements of column in mesh
    elem(idx:(idx+h_elem*2-1),:)=[first'  second'  third']; 
    idx=idx+h_elem*2;
end

%% edge - link between edges and nodes
me_hor=w_elem*(h_elem+1);
me_ver=(w_elem+1)*h_elem;
me_diag=w_elem*h_elem;
me=me_hor+me_ver+me_diag;

tmp=repmat((1:(h_elem+1):(h_elem+1)*(w_elem))',1,h_elem+1)+repmat(0:h_elem,w_elem,1);
edge_hor=[tmp(:) tmp(:)+h_elem+1];
tmp=repmat((1:h_elem)',1,w_elem+1)+repmat((0:(h_elem+1):(h_elem+1)*(w_elem)),h_elem,1);
edge_ver=[tmp(:) tmp(:)+1];
tmp=repmat((1:h_elem)',1,w_elem)+repmat((0:(h_elem+1):(h_elem+1)*(w_elem-1)),h_elem,1);
edge_diag=[tmp(:) tmp(:)+h_elem+2];
edge=[edge_hor;edge_ver;edge_diag];

tri_grid.node=node;
tri_grid.elem=elem;
tri_grid.edge=edge;

%% bdFlag matrix creation
% bdFlag - in each cell (1-bottom,2-right,3-top,4-left)
bdD=cell(4,1);
bdD{1}=1:(h_elem+1):(h_elem+1)*(w_elem+1);
bdD{2}=((h_elem+1)*w_elem+1):(h_elem+1)*(w_elem+1);
bdD{3}=(h_elem+1):(h_elem+1):(h_elem+1)*(w_elem+1);
bdD{4}=1:(h_elem+1);
bdN=cell(4,1);
bdN{1}=1:w_elem;
bdN{2}=me_hor+h_elem*w_elem+(1:h_elem);
bdN{3}=h_elem*w_elem+(1:w_elem);
bdN{4}=me_hor+(1:h_elem);

bound_idx.bdD=bdD;
bound_idx.bdN=bdN;

end