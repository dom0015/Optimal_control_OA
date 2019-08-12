function [N_bound_set] = N_boundary_setting(N_idx,N_val,tri_grid,bound_idx)
if ~iscell(N_val)
    N_val={N_val};
end
n=size(N_idx,1);

tmp_idx=cell(n,1);
tmp_val=cell(n,1);
for i=1:n
    tmp_allidx=bound_idx.bdN{N_idx{i,1}};
    tmp_from=N_idx{i,2}(1);
    tmp_to=N_idx{i,2}(2);
    if N_idx{i,1}==1 || N_idx{i,1}==3
        tmp_e_idx=tri_grid.edge(tmp_allidx,:);
        tmp_n=(tri_grid.node(tmp_e_idx(:,1),:)+tri_grid.node(tmp_e_idx(:,2),:))/2;
        tmp_template=tmp_n(:,1);
    else
        tmp_e_idx=tri_grid.edge(tmp_allidx,:);
        tmp_n=(tri_grid.node(tmp_e_idx(:,1),:)+tri_grid.node(tmp_e_idx(:,2),:))/2;
        tmp_template=tmp_n(:,2);
    end
    
    idx_in=(tmp_template<=tmp_to)&(tmp_template>=tmp_from);
    loc_idx=tmp_allidx(idx_in);
    tmp_edge=tri_grid.edge(loc_idx,:);
    tmp_node=(tri_grid.node(tmp_edge(:,1),:)+tri_grid.node(tmp_edge(:,2),:))/2;
    loc_val=N_val{i}(tmp_node);
    tmp_idx{i}=loc_idx;
    tmp_val{i}=loc_val;
end
idx=cell2mat(tmp_idx')';
val=cell2mat(tmp_val);
m=size(tri_grid.edge,1);
edge_flg=false(m,1);
edge_val=zeros(m,1);

edge_flg(idx)=true;
edge_val(idx)=val;

N_bound_set.edge_flg=edge_flg;
N_bound_set.edge_val=edge_val;
end