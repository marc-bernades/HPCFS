%% Postprocess 3D TGV Inviscid Invariants

% Load results
% From ECCOMAS
KGP      = readtable('output_data_3D_TGV_32x32x32_KGP_125.6637_Time.csv');
KGP_P    = readtable('output_data_3D_TGV_32x32x32_KGP_P_125.6637_Time.csv');
Div      = readtable('output_data_3D_TGV_32x32x32_Div_4.7124_Time.csv');
Div_G    = readtable('output_data_3D_TGV_32x32x32_Div_Gaus_125.6637_Time.csv');
Div_I    = readtable('output_data_3D_TGV_32x32x32_Div_Imp_125.6637_Time.csv');
HLLC     = readtable('output_data_3D_TGV_32x32x32_HLLC_125.6637_Time.csv');
Shima    = readtable('output_data_3D_TGV_32x32x32_Shima_125.6637_Time.csv');

% Normalize time
t_c = 2*pi;

%% Ke evolution
% Main comparison
figure
title('3D TGV Inviscid 32x32x32 @ 20 x t_c')
plot(KGP.t_vec/t_c,KGP.ke_total/KGP.ke_total(1),'LineWidth',2,'LineStyle','-')
hold on
plot(HLLC.t_vec/t_c,HLLC.ke_total/HLLC.ke_total(1),'LineWidth',2,'LineStyle','--')
plot(Div.t_vec/t_c,Div.ke_total/Div.ke_total(1),'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
% plot(DivFiltGaus.t_vec/KGP.t_vec(end),DivFiltGaus.ke_total/DivFiltGaus.ke_total(1),'LineWidth',2,'LineStyle','-.')
plot(Shima.t_vec/t_c,Shima.ke_total/Shima.ke_total(1),'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
% pbaspect([1.8 1 1])
ylim([0.9 1.1])
legend('KGP','HLLC','Div','Shima')
legend('Location','best','box','off')
grid on
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('${{k_e}/{{k_e}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
saveas(gca,'Results/3D_TGV/3D_TGV_Inviscid_ke_Dissipative_20tc','png')

% Dissipative comparison comparison
figure
plot(KGP.t_vec/t_c,KGP.ke_total/KGP.ke_total(1),'LineWidth',2,'LineStyle','-')
hold on
plot(HLLC.t_vec/t_c,HLLC.ke_total/HLLC.ke_total(1),'LineWidth',2,'LineStyle','--','color',[0.4660, 0.6740, 0.1880])
plot(Div_G.t_vec/t_c,Div_G.ke_total/Div_G.ke_total(1),'LineWidth',2,'LineStyle',':')
plot(Div_I.t_vec/t_c,Div_I.ke_total/Div_I.ke_total(1),'LineWidth',2,'LineStyle','-.','color',[0.3010 0.7450 0.9330])
ylim([0.0 1.1])
xlim([0 0.1])
legend('KGP','HLLC','Div + Gaus','Div + Imp')
legend('Location','best','box','off')
grid on
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('${{k_e}/{{k_e}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
saveas(gca,'Results/3D_TGV/3D_TGV_Inviscid_ke_Dissipative_20tc','png')

%% INVARIANTS
% KGP vs KGP w/P
plot_3D_TGV_Invariants_KGP(KGP_P.t_vec, KGP_P.rho_norm, KGP_P.rhou_bar, KGP_P.rhoE_norm, ...
    KGP_P.rhoe_norm, KGP_P.rhosPR_norm, KGP_P.rhoui2_norm, ...
    KGP_P.rhoE2_norm, KGP_P.rhoe2_norm, KGP_P.rhos2_norm,...
    KGP.t_vec, KGP.rho_norm, KGP.rhou_bar, KGP.rhoE_norm, ...
    KGP.rhoe_norm, KGP.rhosPR_norm, KGP.rhoui2_norm, ...
    KGP.rhoE2_norm, KGP.rhoe2_norm, KGP.rhos2_norm);

plot_3D_TGV_Invariants_KGP(KGP.t_vec, KGP.rho_norm, KGP.rhou_bar, KGP.rhoE_norm, ...
    KGP.rhoe_norm, KGP.rhosPR_norm, KGP.rhoui2_norm, ...
    KGP.rhoE2_norm, KGP.rhoe2_norm, KGP.rhos2_norm,...
    Shima.t_vec, Shima.rho_norm, Shima.rhou_bar, Shima.rhoE_norm, ...
    Shima.rhoe_norm, Shima.rhosPR_norm, Shima.rhoui2_norm, ...
    Shima.rhoE2_norm, Shima.rhoe2_norm, Shima.rhos2_norm);
