

% neuman_norm=[];
% dirichlet_norm=[];
% min_eig=[];
% eig_comp_real=[];
% eig_comp_imag=[];
% kkk=0;
% iter_count=[];
% res_res=[];


figure; plot(10.^(-(1:9)),neuman_norm,'-*','LineWidth',3)
grid on
box on
set(gca,'yScale','log')
set(gca,'xScale','log')


title('neuman_norm')
xlabel('value of \beta')
ylabel('neuman_norm')

