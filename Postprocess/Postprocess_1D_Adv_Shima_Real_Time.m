%% Postprocess 3D TGV Inviscid

%% 1D Advective test
% Load results
KGP     = readtable('output_data_1D_Adv_Shima_128x1x1_KGP_0.08_Time.csv');
HLLC_DF      = load('DoubleFlux_128x1x1_KGP.mat');

% Plots > Compare schemes
figure
title('1D Adv N = 128')
Cond = KGP.ke_total>0;
plot(KGP.t_vec(Cond),KGP.ke_total(Cond)/KGP.ke_total(1),'LineWidth',2,'LineStyle',':')
hold on
plot(HLLC_DF.t,HLLC_DF.ke/HLLC_DF.ke(1),'LineWidth',2,'LineStyle',':')
legend('KGP','DF HLLC')
legend('Location','best','box','off')
grid on
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('${{k_e}/{{k_e}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
ylim([0 1.2])
xlim([0 0.08])
pbaspect([1.8 1 1])
saveas(gca,'1D_Adv_Shima_Time','epsc')
