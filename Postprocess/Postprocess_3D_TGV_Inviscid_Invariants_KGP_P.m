%% Postprocess 3D TGV Inviscid Invariants

%% KGP WITH PRESSURE NOVEL SCHEME VS KGP
KGP_P    = readtable('output_data_3D_TGV_32x32x32_KGP_P_125.6637_Time.csv');
KGP      = readtable('output_data_3D_TGV_32x32x32_KGP_125.6637_Time.csv');
% KGP      = readtable('output_data_3D_TGV_32x32x32_KGP_Pt_M_125.6637_Time.csv');
% KGP_P    = readtable('output_data_3D_TGV_32x32x32_KGP_Pt_NCForm_125.6637_Time.csv');
% KGP_P      = readtable('output_data_2D_MixingLayer_Periodic_20x20x1_KGP_P_NCForm_8_Time.csv');
% KGP        = readtable('output_data_2D_MixingLayer_Periodic_20x20x1_KGP_P_GeneralForm_8_Time.csv');
% KGP        = readtable('output_data_2D_MixingLayer_Periodic_20x20x1_KGP_8_Time.csv');
% KGP        = readtable('output_data_2D_MixingLayer_Periodic_Real_20x20x1_KGP_P_1.2566_Time.csv');

% KGP        = readtable('output_data_2D_Vortex_40x40x1_KGP-Et_4_Time.csv');
% KGP_P        = readtable('output_data_2D_Vortex_40x40x1_KGP-Pt_4_Time.csv');

% Plot features
c_label  = {'KGP-Pt','KGP-Et'}; % Legend
% c_label  = {'NC-2','RG'}; % Legend

t_c      = 2*pi; % Characteristic time
t_end    = 20.0; % Number of cycles to plot


% Plots > Compare schemes 20FTT
figure
plot(KGP_P.t_vec/t_c,KGP_P.ke_total/KGP_P.ke_total(1),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(KGP.t_vec/t_c,KGP.ke_total/KGP.ke_total(1),'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
pbaspect([1.8 1 1])
legend(c_label)
legend('Location','best','box','off')
grid on
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('${{k_e}/{{k_e}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
xlim([0 t_end])
saveas(gca,'Results/3D_TGV/3D_TGV_Inviscid_ke_KGP_wP','png')


plot_3D_TGV_Invariants_KGP(KGP_P.t_vec, KGP_P.rho_norm, KGP_P.rhou_bar, KGP_P.rhoE_norm, ...
    KGP_P.P_norm, KGP_P.rhosPR_norm, KGP_P.rhoui2_norm, ...
    KGP_P.rhoE2_norm, KGP_P.rhoe2_norm, KGP_P.rhosPR2_norm,...
    KGP.t_vec, KGP.rho_norm, KGP.rhou_bar, KGP.rhoE_norm, ...
    KGP.P_norm, KGP.rhosPR_norm, KGP.rhoui2_norm, ...
    KGP.rhoE2_norm, KGP.rhoe2_norm, KGP.rhosPR2_norm,t_c,c_label,t_end);


%% KGP VALIDATION 32^3 vs 64^3
%% KGP WITH PRESSURE NOVEL SCHEME VS KGP
KGP_P_64      = readtable('output_data_3D_TGV_16x16x16_KGP_P_125.6637_Time.csv');

% Plot features
c_label  = {'KGP w/P 32^3','KGP w/P 16^3'}; % Legend
t_c      = 2*pi; % Characteristic time
t_end    = 20; % Number of cycles to plot


% Plots > Compare schemes 20FTT
figure
plot(KGP_P.t_vec/t_c,KGP_P.ke_total/KGP_P.ke_total(1),'LineWidth',2,'LineStyle','-')
hold on
plot(KGP_P_64.t_vec/t_c,KGP_P_64.ke_total/KGP_P_64.ke_total(1),'LineWidth',2,'LineStyle','--')
pbaspect([1.8 1 1])
legend(c_label)
legend('Location','best','box','off')
grid on
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('${{k_e}/{{k_e}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlim([0 t_end])
saveas(gca,'Results/3D_TGV/3D_TGV_Inviscid_ke_KGP_wP_convergence','png')


plot_3D_TGV_Invariants_KGP(KGP_P.t_vec, KGP_P.rho_norm, KGP_P.rhou_bar, KGP_P.rhoE_norm, ...
    KGP_P.P_norm, KGP_P.rhosPR_norm, KGP_P.rhoui2_norm, ...
    KGP_P.rhoE2_norm, KGP_P.rhoe2_norm, KGP_P.rhosPR2_norm,...
    KGP_P_64.t_vec, KGP_P_64.rho_norm, KGP_P_64.rhou_bar, KGP_P_64.rhoE_norm, ...
    KGP_P_64.P_norm, KGP_P_64.rhosPR_norm, KGP_P_64.rhoui2_norm, ...
    KGP_P_64.rhoE2_norm, KGP_P_64.rhoe2_norm, KGP_P_64.rhosPR2_norm,t_c,c_label,t_end);

%% FILTER COMPARISON
%% 32^3 CFL 0.3
Div      = readtable('output_data_3D_TGV_32x32x32_Div_4.7124_Time.csv');
Div_F2   = readtable('output_data_3D_TGV_32x32x32_Div_F2_125.6637_Time.csv');
Div_F4   = readtable('output_data_3D_TGV_32x32x32_Div_F4_125.6637_Time.csv');

% Div blows up at ~3s. F4 follows div up to 3s.
% Alpha values of 0.495. Alpha 0.499 less dissipation for F2, but F4 more
% stable to garantee margin to 0.5 limit

% Plot features
c_label  = {'Div','Div w/F2','Div w/F4'}; % Legend
t_c      = 2*pi; % Characteristic time
t_end    = 0.5;  % Number of cycles to plot


% Plots > Compare schemes 20FTT
figure
plot(Div.t_vec/t_c,Div.ke_total/Div.ke_total(1),'LineWidth',2,'LineStyle','-')
hold on
plot(Div_F2.t_vec/t_c,Div_F2.ke_total/Div_F2.ke_total(1),'LineWidth',2,'LineStyle','--')
plot(Div_F4.t_vec/t_c,Div_F4.ke_total/Div_F4.ke_total(1),'LineWidth',2,'LineStyle','-.')
pbaspect([1.8 1 1])
legend(c_label)
legend('Location','best','box','off')
grid on
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('${{k_e}/{{k_e}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlim([0 t_end])
saveas(gca,'Results/3D_TGV/3D_TGV_Inviscid_ke_Filter','png')

%% FILTER COMPARISON
%% 64^3 CFL 0.3
Div      = readtable('output_data_3D_TGV_64x64x64_Div_3.1416_Time.csv');
Div_F2   = readtable('output_data_3D_TGV_64x64x64_Div_F2_6.2832_Time.csv');
Div_F4   = readtable('output_data_3D_TGV_64x64x64_Div_F4_6.2832_Time.csv');

% Plot features
c_label  = {'Div','Div w/F2','Div w/F4'}; % Legend
t_c      = 2*pi; % Characteristic time
t_end    = 0.5; % Number of cycles to plot


% Plots > Compare schemes 20FTT
figure
plot(Div.t_vec/t_c,Div.ke_total/Div.ke_total(1),'LineWidth',2,'LineStyle','-')
hold on
plot(Div_F2.t_vec/t_c,Div_F2.ke_total/Div_F2.ke_total(1),'LineWidth',2,'LineStyle','--')
plot(Div_F4.t_vec/t_c,Div_F4.ke_total/Div_F4.ke_total(1),'LineWidth',2,'LineStyle','-.')
pbaspect([1.8 1 1])
legend(c_label)
legend('Location','best','box','off')
grid on
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('${{k_e}/{{k_e}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlim([0 t_end])
saveas(gca,'Results/3D_TGV/3D_TGV_Inviscid_ke_Filter','png')

%% KGP WITH PRESSURE OTHER SCHEMES
DATA       = readtable('output_data_3D_TGV_32x32x32_KGP_125.6637_Time.csv');
% DATA2      = readtable('output_data_3D_TGV_32x32x32_Div_4.7124_Time.csv');
DATA2      = readtable('output_data_3D_TGV_32x32x32_Div_F4_125.6637_Time.csv');

% Plot features
c_label  = {'KGP','Div wF4'}; % Legend
t_c      = 2*pi; % Characteristic time
t_end    = 1; % Number of cycles to plot


plot_3D_TGV_Invariants_KGP(DATA.t_vec, DATA.rho_norm, DATA.rhou_bar, DATA.rhoE_norm, ...
    DATA.P_norm, DATA.rhosPR_norm, DATA.rhoui2_norm, ...
    DATA.rhoE2_norm, DATA.rhoe2_norm, DATA.rhosPR2_norm,...
    DATA2.t_vec, DATA2.rho_norm, DATA2.rhou_bar, DATA2.rhoE_norm, ...
    DATA2.P_norm, DATA2.rhosPR_norm, DATA2.rhoui2_norm, ...
    DATA2.rhoE2_norm, DATA2.rhoe2_norm, DATA2.rhosPR2_norm,t_c,c_label,t_end);