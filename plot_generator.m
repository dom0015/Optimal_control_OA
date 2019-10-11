

% neuman_norm=[];
% dirichlet_norm=[];
% min_eig=[];
% eig_comp_real=[];
% eig_comp_imag=[];
% kkk=0;
% iter_count=[];
% res_res=[];


figure; plot(10.^(-(1:9)),dirichlet_norm,'-*','LineWidth',3)
grid on
box on
set(gca,'yScale','log')
set(gca,'xScale','log')


title('dirichlet_norm')
xlabel('value of \beta')
ylabel('dirichlet_norm')

