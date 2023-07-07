%% Postprocess 1D Adv Inviscid Shima
close all

%% 1D Advective test
%% Load results
KGP        = readtable('output_data_1D_Adv_Shima_41x1x1_KGP_0.01.csv');
Div        = readtable('output_data_1D_Adv_Shima_41x1x1_Div_0.01.csv');

% DivFilt_G  = readtable('output_data_1D_Adv_Shima_41x1x1_Div_Gau_0.0033243.csv');
DivFilt_I  = readtable('output_data_1D_Adv_Shima_41x1x1_Div_Imp_0.01.csv');
% DivFilt_I2  = readtable('output_data_1D_Adv_Shima_41x1x1_Div_Filt_Imp_0.01.csv');

HLLC       = readtable('output_data_1D_Adv_Shima_41x1x1_HLLC_0.01.csv');
Shima      = readtable('output_data_1D_Adv_Shima_41x1x1_Shima_0.01.csv');
KGP_P      = readtable('output_data_1D_Adv_Shima_41x1x1_KGP_P_0.01.csv');
% PEP_RGo2   = readtable('output_data_1D_Adv_Shima_41x1x1_Real_PEP-RG_0.01.csv');
KGP_Et_CR    = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-Et-CR_0.01.csv');
KGP_et_CR    = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-et-CR_0.01.csv');

%% Final implementation results
Div      = readtable('output_data_1D_Adv_Shima_41x1x1_Real_Div_0.01.csv');
Div_F4   = readtable('output_data_1D_Adv_Shima_41x1x1_Real_Div_F4_0.01.csv');
Div_F2   = readtable('output_data_1D_Adv_Shima_41x1x1_Real_Div_F2_0.01.csv');
Hybrid   = readtable('output_data_1D_Adv_Shima_41x1x1_Real_Hybrid_0.01.csv');
WENO5    = readtable('output_data_1D_Adv_Shima_41x1x1_Real_WENO5_1.csv');
DivPt_F4 = readtable('output_data_1D_Adv_Shima_41x1x1_Real_D-Pt_F4_0.01.csv');

% Index only unique Y and Z
Index_Z = unique(KGP.Z); Index_Z = find(KGP.Z == Index_Z(2));
Index   = Index_Z(5:3:end-3);

% Normalize
P_0   = 5E6;
u_0   = 1;
rho_0 = 95.496;

% Reference solution
rho_min     = 56.9;                   % Min density non-linear
rho_max     = 793.1;                  % Max density non-linear
rho_ref = (rho_min + rho_max)/2  + (rho_max - rho_min)/2*sin(2*pi*KGP.X(Index));

% Data trim for marker dashed line
HLLC_X   = HLLC.X(Index);
HLLC_u   = HLLC.u(Index);
HLLC_P   = HLLC.P(Index);
HLLC_rho = HLLC.rho(Index);

% DF load
DF_HLLC = load('DoubleFlux_41x1x1_HLLC.mat');
DF_KGP = load('DoubleFlux_41x1x1_KGP.mat');

%% 'KGP w/P','Shima','DF w/HLLC','DF w/KGP'
% Plots > Compare schemes
figure
% subplot(1,3,1)
% title('1D Inviscid Advective N = 41')
% plot(KGP_P.X(Index),KGP_P.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-')
plot(KGP.X(Index),rho_ref/rho_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
plot(KGP_P.X(Index),KGP_P.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(KGP.X(Index),KGP.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(KGP_Et_CR.X(Index),KGP_Et_CR.rho(Index)/rho_0,'LineWidth',2, 'LineStyle',':','color',[0.3010 0.7450 0.9330])
plot(KGP_et_CR.X(Index),KGP_et_CR.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-.','color',[0, 0, 1])
% plot(HLLC_X(1:2:end),HLLC_rho(1:2:end)/rho_0,'--o','LineWidth',2,'MarkerSize',3)
% plot(Div.X(Index),Div.rho(Index)/rho_0,'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
% plot(DivFilt_I.X(Index),DivFilt_I.rho(Index)/rho_0,'LineWidth',2,'LineStyle','-.')
plot(Shima.X(Index),Shima.rho(Index)/rho_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660 0.6740 0.1880])
plot(DF_KGP.x, DF_KGP.V(1,:)/rho_0,'LineWidth',2,'LineStyle','-.','color',[0.8500 0.3250 0.0980])
plot(DF_HLLC.x, DF_HLLC.V(1,:)/rho_0,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
plot(Div_F4.X(Index),Div_F4.rho(Index)/rho_0,'LineWidth',2,'LineStyle','--','Color',[0.4940 0.1840 0.5560])
plot(DivPt_F4.X(Index(1:2:end)),DivPt_F4.rho(Index(1:2:end))/rho_0,'o','LineWidth',1.5,'MarkerSize',8,'color',[0.9290 0.6940 0.1250])
% plot(PEP_RGo2.X(Index),PEP_RGo2.rho(Index)/rho_0,'x','LineWidth',2,'MarkerSize',6,'color','black')

% plot(Hybrid.X(Index),Hybrid.rho(Index)/rho_0,'LineWidth',2,'marker','o','color',[0.8500 0.3250 0.0980])
legend('Exact','KGP-Pt','KGP-Et','KGP-Et-CR','KGP-et-CR','PEP-IG','KGP-Df','UB-Df','D+F4','D-Pt+F4')

% legend('KGP-Pt','Shima','KGP-DF','HLLC-DF')
legend('Location','best','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{\rho}/{{\rho}_{0}}}$','interpreter','latex')
ylim([0 9])
% ylabel('\rho [kg/m^3]')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

% Inset
axes('position',[.205 .2 .35 .3])
plot(KGP.X(Index),rho_ref*0 + P_0/P_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
% plot(KGP_P.X(Index),KGP_P.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-')
plot(KGP.X(Index),rho_ref/rho_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
plot(KGP_P.X(Index),KGP_P.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(KGP.X(Index),KGP.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(KGP_Et_CR.X(Index),KGP_Et_CR.rho(Index)/rho_0,'LineWidth',2, 'LineStyle',':','color',[0.3010 0.7450 0.9330])
plot(KGP_et_CR.X(Index),KGP_et_CR.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-.','color',[0, 0, 1])
% plot(HLLC_X(1:2:end),HLLC_rho(1:2:end)/rho_0,'--o','LineWidth',2,'MarkerSize',3)
% plot(Div.X(Index),Div.rho(Index)/rho_0,'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
% plot(DivFilt_I.X(Index),DivFilt_I.rho(Index)/rho_0,'LineWidth',2,'LineStyle','-.')
plot(Shima.X(Index),Shima.rho(Index)/rho_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660 0.6740 0.1880])
plot(DF_KGP.x, DF_KGP.V(1,:)/rho_0,'LineWidth',2,'LineStyle','-.','color',[0.8500 0.3250 0.0980])
plot(DF_HLLC.x, DF_HLLC.V(1,:)/rho_0,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
plot(Div_F4.X(Index),Div_F4.rho(Index)/rho_0,'LineWidth',2,'LineStyle','--','Color',[0.4940 0.1840 0.5560])
plot(DivPt_F4.X(Index(1:2:end)),DivPt_F4.rho(Index(1:2:end))/rho_0,'o','LineWidth',1.5,'MarkerSize',8,'color',[0.9290 0.6940 0.1250])

xlim([0.2 0.3])
ylim([max(rho_ref/rho_0)*0.98 max(rho_ref/rho_0)*1.006])

set(gca,'linewidth',1)
set(gca,'fontsize',10)

saveas(gca,'Results/1D_Adv_Shima/Real/41x1x1/1D_Shima_Real_rho','epsc')

figure
% subplot(1,3,2)
% plot(KGP.X(Index),KGP.P(Index)/P_0,'LineWidth',1.5, 'LineStyle','-')
% subplot(1,3,2)
% plot(KGP.X(Index),KGP.P(Index)/P_0,'LineWidth',1.5, 'LineStyle','-')
plot(KGP.X(Index),rho_ref*0 + P_0/P_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
plot(KGP_P.X(Index),KGP_P.P(Index)/P_0,'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(KGP.X(Index),KGP.P(Index)/P_0,'LineWidth',2, 'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(KGP_Et_CR.X(Index),KGP_Et_CR.P(Index)/P_0,'LineWidth',2, 'LineStyle',':','color',[0.3010 0.7450 0.9330])
plot(KGP_et_CR.X(Index),KGP_et_CR.P(Index)/P_0,'LineWidth',2, 'LineStyle','-.','color',[0, 0, 1])
plot(Shima.X(Index),Shima.P(Index)/P_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660 0.6740 0.1880])
plot(DF_KGP.x, DF_KGP.V(3,:)/P_0,'LineWidth',2,'LineStyle','-.','color',[0.8500 0.3250 0.0980])
plot(DF_HLLC.x, DF_HLLC.V(3,:)/P_0,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
plot(Div_F4.X(Index),Div_F4.P(Index)/P_0,'LineWidth',2,'LineStyle','--','color',[0.4940 0.1840 0.5560])
plot(DivPt_F4.X(Index(1:2:end)),DivPt_F4.P(Index(1:2:end))/P_0,'o','LineWidth',1.5,'MarkerSize',8,'color',[0.9290 0.6940 0.1250])
% plot(PEP_RGo2.X(Index),PEP_RGo2.P(Index)/P_0,'x','LineWidth',2,'MarkerSize',6,'color',[0 0.4470 0.7410])
% plot(KGP_snew.X(Index),KGP_new.P(Index)/P_0,'LineStyle','-.','LineWidth',2,'MarkerSize',6,'color',[0 0.4470 0.7410])
% legend('Exact','KGP-Pt','KGP-Et','KGP-Et-CR','KGP-et-CR','PEP-IG','KGP-Df','UB-Df','D+F4','D-Pt+F4')

ylim([0.992,1.008])
yticklabels(0.992:0.002:1.008)

% legend('Exact','KGP-Pt','Shima','KGP-Df','UB-Df','D+F4')
% legend('Location','northwest','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{P}/{{P}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

axes('position',[.275 .2 .35 .3])
plot(KGP.X(Index),rho_ref*0 + P_0/P_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
plot(KGP_P.X(Index),KGP_P.P(Index)/P_0,'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(KGP_Et_CR.X(Index),KGP_Et_CR.P(Index)/P_0,'LineWidth',2, 'LineStyle',':','color',[0.3010 0.7450 0.9330])
plot(KGP_et_CR.X(Index),KGP_et_CR.P(Index)/P_0,'LineWidth',2, 'LineStyle','-.','color',[0, 0, 1])
plot(Shima.X(Index),Shima.P(Index)/P_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660 0.6740 0.1880])
plot(DF_KGP.x, DF_KGP.V(3,:)/P_0,'LineWidth',2,'LineStyle','-.','color',[0.8500 0.3250 0.0980])
plot(DF_HLLC.x, DF_HLLC.V(3,:)/P_0,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
plot(DivPt_F4.X(Index(1:4:end)),DivPt_F4.P(Index(1:4:end))/P_0,'o','LineWidth',1.5,'MarkerSize',8,'color',[0.9290 0.6940 0.1250])
% plot(PEP_RGo2.X(Index),PEP_RGo2.P(Index)/P_0,'x','LineWidth',2,'MarkerSize',6,'color',[0 0.4470 0.7410])
% plot(KGP_new.X(Index),KGP_new.P(Index)/P_0,'LineStyle','-.','LineWidth',2,'MarkerSize',6,'color',[0 0.4470 0.7410])

ylim([0.9998,1.0002])
set(gca,'linewidth',1)
set(gca,'fontsize',10)

saveas(gca,'Results/1D_Adv_Shima/Real/41x1x1/1D_Shima_Real_P','epsc')

figure
% subplot(1,3,3)
plot(KGP.X(Index),rho_ref*0 + u_0/u_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
plot(KGP_P.X(Index),KGP_P.u(Index)/u_0,'LineWidth',2, 'LineStyle','-','color',[0 0.4470 0.7410])
plot(KGP.X(Index),KGP.u(Index)/u_0,'LineWidth',2, 'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(KGP_Et_CR.X(Index),KGP_Et_CR.u(Index)/u_0,'LineWidth',2, 'LineStyle',':','color',[0.3010 0.7450 0.9330])
plot(KGP_et_CR.X(Index),KGP_et_CR.u(Index)/u_0,'LineWidth',2, 'LineStyle','-.','color',[0, 0, 1])
plot(Shima.X(Index),Shima.u(Index)/u_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660 0.6740 0.1880])
plot(DF_KGP.x, DF_KGP.V(2,:)/u_0,'LineWidth',2,'LineStyle','-.','color',[0.8500 0.3250 0.0980])
plot(DF_HLLC.x, DF_HLLC.V(2,:)/u_0,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
plot(Div_F4.X(Index),Div_F4.u(Index)/u_0,'LineWidth',2,'LineStyle','--','color',[0.4940 0.1840 0.5560])


legend('Exact','KGP-Pt','KGP-Et','KGP-Et-CR','KGP-et-CR','PEP-IG','KGP-Df','UB-Df','D+F4','D-Pt+F4')
legend('Location','southwest','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{u}/{{u}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,'Results/1D_Adv_Shima/Real/41x1x1/1D_Shima_Real_u','epsc')

%% 'KGP w/P','DF w/HLLC','Hybrid','Div w/F4'
% Plots > Compare schemes
figure
% subplot(1,3,1)
% title('1D Inviscid Advective N = 41')
% plot(KGP_P.X(Index),KGP_P.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-')
hold on
plot(KGP_P.X(Index),KGP_P.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-')
plot(DF_HLLC.x, DF_HLLC.V(1,:)/rho_0,'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
plot(Hybrid.X(Index),Hybrid.rho(Index)/rho_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Div_F4.X(Index),Div_F4.rho(Index)/rho_0,'LineWidth',2,'LineStyle','-.','color',[0.6350 0.0780 0.1840])

legend('KGP w/P','DF w/HLLC','Hybrid','Div w/F4')
legend('Location','best','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{\rho}/{{\rho}_{0}}}$','interpreter','latex')
% ylabel('\rho [kg/m^3]')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
saveas(gca,'Results/1D_Adv_Shima/Real/41x1x1/1D_Shima_Real_rho_2','png')

figure
% subplot(1,3,2)
hold on
plot(KGP_P.X(Index),KGP_P.P(Index)/P_0,'LineWidth',2, 'LineStyle','-')
plot(DF_HLLC.x, DF_HLLC.V(3,:)/P_0,'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
plot(Hybrid.X(Index),Hybrid.P(Index)/P_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Div_F4.X(Index),Div_F4.P(Index)/P_0,'LineWidth',2,'LineStyle','-.','color',[0.6350 0.0780 0.1840])

legend('KGP w/P','DF w/HLLC','Hybrid','Div w/F4')
legend('Location','best','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{P}/{{P}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
saveas(gca,'Results/1D_Adv_Shima/Real/41x1x1/1D_Shima_Real_P_2','png')

figure
% subplot(1,3,3)
hold on
plot(KGP_P.X(Index),KGP_P.u(Index)/u_0,'LineWidth',2, 'LineStyle','-')
plot(DF_HLLC.x, DF_HLLC.V(2,:)/u_0,'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
plot(Hybrid.X(Index),Hybrid.u(Index)/u_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Div_F4.X(Index),Div_F4.u(Index)/u_0,'LineWidth',2,'LineStyle','-.','color',[0.6350 0.0780 0.1840])

legend('KGP w/P','DF w/HLLC','Hybrid','Div w/F4')

legend('Location','southwest','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{u}/{{u}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
saveas(gca,'Results/1D_Adv_Shima/Real/41x1x1/1D_Shima_Real_u_2','png')