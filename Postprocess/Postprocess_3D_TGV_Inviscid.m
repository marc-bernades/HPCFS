%% Postprocess 3D TGV Inviscid

% Load results
% From ECCOMAS
KGP     = readtable('output_data_3D_TGV_32x32x32_KGP_125.6637_Time.csv');
HLLC    = readtable('output_data_3D_TGV_32x32x32_HLLC_Time.csv');
Div     = readtable('output_data_3D_TGV_32x32x32_Div_5-3_Time.csv');
DivFilt = readtable('output_data_3D_TGV_32x32x32_DivFilt_5-1_Time.csv');
DivFiltGaus = readtable('output_data_3D_TGV_32x32x32_DivFiltGaus_Time.csv');
DivFiltImp  = readtable('output_data_3D_TGV_32x32x32_Div_Imp_31.4159_Time.csv');
% CTR
Shima   = readtable('output_data_3D_TGV_32x32x32_Shima_Time.csv');

% Reference time
t_c = 2*pi;

% Plots > Compare schemes
figure
title('3D TGV Inviscid 32x32x32 @ 10s')
plot(KGP.t_vec/t_c,KGP.ke_total/KGP.ke_total(1),'LineWidth',2,'LineStyle','-')
hold on
plot(HLLC.t_vec(10:50:end)/t_c,HLLC.ke_total(10:50:end)/HLLC.ke_total(1),'--o','LineWidth',2,'MarkerSize',3)
plot(Div.t_vec/t_c,Div.ke_total/Div.ke_total(1),'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
plot(DivFiltGaus.t_vec/t_c,DivFiltGaus.ke_total/DivFiltGaus.ke_total(1),'LineWidth',2,'LineStyle','-.')
plot(Shima.t_vec/t_c,Shima.ke_total/Shima.ke_total(1),'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])

legend('KGP','HLLC','Div','Div + Filter','Shima')
legend('Location','best','box','off')
grid on
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('${{k_e}/{{k_e}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
ylim([0 1.2])
pbaspect([1.8 1 1])
saveas(gca,'Results/3D_TGV/3D_TGV_Inviscid','epsc')
