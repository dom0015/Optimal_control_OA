function [] = plot_grid(tri_grid)
%PLOT_GRID Summary of this function goes here
%   Detailed explanation goes here
node=tri_grid.node;
elem=tri_grid.elem;
edge=tri_grid.edge;
txt_size=10;
figure
hold on
plot(node(:,1),node(:,2),'g.','MarkerSize',20);
for i=1:size(node,1)
    text(node(i,1),node(i,2),num2str(i),'Color','black','FontSize',txt_size);
end
for i=1:size(elem,1)
    tmp_p=node(elem(i,:),:);
    x=sum(tmp_p(:,1))/3;
    y=sum(tmp_p(:,2))/3;
    text(x,y,num2str(i),'Color','red','FontSize',txt_size);
end
for i=1:size(edge,1)
    tmp_p=node(edge(i,:),:);
    plot(tmp_p(:,1),tmp_p(:,2),'k-');
    x=sum(tmp_p(:,1))/2;
    y=sum(tmp_p(:,2))/2;
    text(x,y,num2str(i),'Color','blue','FontSize',txt_size);
end
end

