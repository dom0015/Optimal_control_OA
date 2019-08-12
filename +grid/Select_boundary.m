function [Boundary_edges] = Select_boundary(B_idx,tri_grid,bound_idx)
n=size(B_idx,1);

tmp_idx=cell(n,1);
for i=1:n
    tmp_allidx=bound_idx.bdN{B_idx{i,1}};
    tmp_from=B_idx{i,2}(1);
    tmp_to=B_idx{i,2}(2);
    %tmp_template=linspace(0+0.5/length(tmp_allidx),1-0.5/length(tmp_allidx),length(tmp_allidx));
    if B_idx{i,1}==1 || B_idx{i,1}==3
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
    tmp_idx{i}=loc_idx;
end
idx=cell2mat(tmp_idx')';
m=size(tri_grid.edge,1);
edge_flg=false(m,1);

edge_flg(idx)=true;

Boundary_edges=edge_flg;
end