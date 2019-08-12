function [D_bound_set] = D_boundary_setting(D_idx,D_val,tri_grid,bound_idx)
if ~iscell(D_val)
    D_val={D_val};
end
n=size(D_idx,1);

tmp_idx=cell(n,1);
tmp_val=cell(n,1);
for i=1:n
    tmp_allidx=bound_idx.bdD{D_idx{i,1}};
    tmp_from=D_idx{i,2}(1);
    tmp_to=D_idx{i,2}(2);
    if D_idx{i,1}==1 || D_idx{i,1}==3
        tmp_template=tri_grid.node(tmp_allidx,1);
    else
        tmp_template=tri_grid.node(tmp_allidx,2);
    end
    idx_in=(tmp_template<=tmp_to)&(tmp_template>=tmp_from);
    loc_idx=tmp_allidx(idx_in);
    tmp_node=tri_grid.node(loc_idx,:);
    loc_val=D_val{i}(tmp_node);
    tmp_idx{i}=loc_idx;
    tmp_val{i}=loc_val;
end
idx=cell2mat(tmp_idx')';
val=cell2mat(tmp_val);
m=size(tri_grid.node,1);
node_flg=false(m,1);
node_val=zeros(m,1);

node_flg(idx)=true;
node_val(idx)=val;

D_bound_set.node_flg=node_flg;
D_bound_set.node_val=node_val;
end