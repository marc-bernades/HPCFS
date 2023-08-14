%% Postprocess 1D Adv Inviscid Shima

%% 1D Advective test
% Load results
KGP        = readtable('output_data_1D_Adv_Shima_41x1x1Ideal_KGP.csv');
KGP_P      = readtable('output_data_1D_Adv_Shima_41x1x1_KGP_P_11.csv');
HLLC       = readtable('output_data_1D_Adv_Shima_41x1x1Ideal_HLLC.csv');
Div        = readtable('output_data_1D_Adv_Shima_41x1x1Ideal_Div.csv');
DivFilt_G  = readtable('output_data_1D_Adv_Shima_41x1x1Ideal_Div_Gau.csv');
% DivFilt_I  = readtable('output_data_1D_Adv_Shima_41x1x1Ideal_Div_Imp.csv');
Shima      = readtable('output_data_1D_Adv_Shima_41x1x1Ideal_Shima.csv');

% Load additional cases (CTR)
Div2       = readtable('output_data_1D_Adv_Shima_41x1x1_Div_11.csv');
Div_F2     = readtable('output_data_1D_Adv_Shima_41x1x1_Div_F2_11.csv');
Div_F4     = readtable('output_data_1D_Adv_Shima_41x1x1_Div_F4_11.csv');
% Div_F4     = readtable('output_data_1D_Adv_Shima_41x1x1_D_Pt_F4_11.csv');

Div_G      = readtable('output_data_1D_Adv_Shima_41x1x1_Div_G_11.csv');
Hybrid     = readtable('output_data_1D_Adv_Shima_41x1x1_Hybrid_11.csv');
Hybrid     = readtable('output_data_1D_Adv_Shima_41x1x1_H_11.csv');

Hybrid_HLLC  = readtable('output_data_1D_Adv_Shima_41x1x1_Hybrid_HLLC_11.csv');
WENO5      = readtable('output_data_1D_Adv_Shima_41x1x1_WENO5_11.csv');


% Index only unique Y and Z
Index_Z = unique(KGP.Z); Index_Z = find(KGP.Z == Index_Z(2));
Index   = Index_Z(5:3:end-3);

% Normalize
P_0   = 1;
u_0   = 1;
rho_0 = 2;

% Reference solution
rho_max = 3;
rho_min = 1;
rho_ref = (rho_min + rho_max)/2  + (rho_max - rho_min)/2*sin(2*pi*KGP.X(Index));

% DAta trim for marker dashed line
HLLC_X   = HLLC.X(Index);
HLLC_u   = HLLC.u(Index);
HLLC_P   = HLLC.P(Index);
HLLC_rho = HLLC.rho(Index);

WENO5_X   = WENO5.X(Index);
WENO5_u   = WENO5.u(Index);
WENO5_P   = WENO5.P(Index);
WENO5_rho = WENO5.rho(Index);

Hybrid_X   = Hybrid.X(Index);
Hybrid_u   = Hybrid.u(Index);
Hybrid_P   = Hybrid.P(Index);
Hybrid_rho = Hybrid.rho(Index);

%% 'KGP','HLLC','Div','Shima'
% Plots > Compare schemes
figure
% subplot(1,3,1)
title('1D Inviscid Advective N = 41 @ 11xFTT')
plot(KGP.X(Index),rho_ref/rho_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
plot(KGP_P.X(Index),KGP_P.rho(Index)/rho_0','LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(KGP.X(Index),KGP.rho(Index)/rho_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Shima.X(Index),Shima.rho(Index)/rho_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660, 0.6740, 0.1880])
plot(Div.X(Index),Div.rho(Index)/rho_0,'o','LineWidth',2,'MarkerSize',8,'color',[0.4940 0.1840 0.5560])
plot(HLLC.X(Index),HLLC.rho(Index)/rho_0,'LineWidth',2,'LineStyle','-.','color',[0.6350 0.0780 0.1840])
plot(Hybrid.X(Index),Hybrid.rho(Index)/rho_0,'+','LineWidth',2,'MarkerSize',6,'color',[0.8500 0.3250 0.0980])

legend('Exact','KGP-Pt','KGP-Et','PEP-IG','D','UB','H')
legend('Location','northeast','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{\rho}/{{\rho}_{0}}}$','interpreter','latex')
% ylabel('\rho [kg/m^3]')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,'Results/1D_Adv_Shima/Ideal/1D_Shima_Ideal_rho','epsc')

figure
% subplot(1,3,2)
plot(KGP.X(Index),rho_ref*0+P_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
plot(KGP_P.X(Index),KGP_P.P(Index)/P_0','LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(KGP.X(Index),KGP.P(Index)/P_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Shima.X(Index),Shima.P(Index)/P_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660, 0.6740, 0.1880])
plot(Div.X(Index),Div.P(Index)/P_0,'o','LineWidth',2,'MarkerSize',8,'color',[0.4940 0.1840 0.5560])
plot(HLLC.X(Index),HLLC.P(Index)/P_0,'LineWidth',2,'LineStyle','-.','color',[0.6350 0.0780 0.1840])
% plot(Hybrid.X(Index),Hybrid.P(Index)/P_0,'+','LineWidth',2,'MarkerSize',6,'color',[0.8500 0.3250 0.0980])

% legend('Exact','KGP-Pt','KGP-Et','Shima','D','UB','H')
% legend('Location','best','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{P}/{{P}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

axes('position',[.50 .2 .35 .3])
plot(KGP.X(Index),rho_ref*0+P_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
% plot(KGP_P.X(Index),KGP_P.P(Index)/P_0','LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
% plot(KGP.X(Index),KGP.P(Index)/P_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
% plot(Shima.X(Index),Shima.P(Index)/P_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660, 0.6740, 0.1880])
% plot(Div.X(Index),Div.P(Index)/P_0,'o','LineWidth',2,'MarkerSize',8,'color',[0.4940 0.1840 0.5560])
% plot(HLLC.X(Index),HLLC.P(Index)/P_0,'LineWidth',2,'LineStyle','-.','color',[0.6350 0.0780 0.1840])
plot(Hybrid.X(Index),Hybrid.P(Index)/P_0,'+','LineWidth',2,'MarkerSize',6,'color',[0.8500 0.3250 0.0980])
ylim([0.99997,1.00001])
set(gca,'linewidth',1)
set(gca,'fontsize',10)


saveas(gca,'Results/1D_Adv_Shima/Ideal/1D_Shima_Ideal_P','epsc')

figure
% subplot(1,3,3)
plot(KGP.X(Index),rho_ref*0+u_0,'LineWidth',2, 'LineStyle',':','Color','black')
hold on
plot(KGP_P.X(Index),KGP_P.u(Index)/u_0','LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(KGP.X(Index),KGP.u(Index)/u_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Shima.X(Index),Shima.u(Index)/u_0,'x','LineWidth',2,'MarkerSize',6,'color',[0.4660, 0.6740, 0.1880])
plot(Div.X(Index),Div.u(Index)/u_0,'o','LineWidth',2,'MarkerSize',8,'color',[0.4940 0.1840 0.5560])
plot(HLLC.X(Index),HLLC.u(Index)/u_0,'LineWidth',2,'LineStyle','-.','color',[0.6350 0.0780 0.1840])
plot(Hybrid.X(Index),Hybrid.u(Index)/u_0,'+','LineWidth',2,'MarkerSize',6,'color',[0.8500 0.3250 0.0980])

legend('Exact','KGP-Pt','KGP-Et','PEP-IG','D','UB','H')
legend('Location','southwest','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{u}/{{u}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,'Results/1D_Adv_Shima/Ideal/1D_Shima_Ideal_u','epsc')

%% 'Div','Div wF4','Hybrid','WENO5'
% Plots > Compare schemes
figure
% subplot(1,3,1)
title('1D Inviscid Advective N = 41 @ 11xFTT')
plot(Div2.X(Index),Div2.rho(Index)/rho_0,'LineWidth',2, 'LineStyle','-')
hold on
plot(Div_F4.X(Index),Div_F4.rho(Index)/rho_0,'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
% plot(DivFilt_G.X(Index),DivFilt_G.rho(Index)/rho_0,'LineWidth',2,'LineStyle','-.')
% plot(Div_F2.X(Index),Div_F2.rho(Index)/rho_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Hybrid.X(Index),Hybrid.rho(Index)/rho_0,'LineWidth',2,'LineStyle','-.','color',[0.3010 0.7450 0.9330])
plot(WENO5_X(1:2:end),WENO5_rho(1:2:end)/rho_0,'--o','LineWidth',2,'MarkerSize',3,'color',[0.6350 0.0780 0.1840])

legend('Div','Div w/F4','Hybrid','WENO5')
legend('Location','best','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{\rho}/{{\rho}_{0}}}$','interpreter','latex')
% ylabel('\rho [kg/m^3]')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
saveas(gca,'Results/1D_Adv_Shima/Ideal/1D_Shima_Ideal_rho_2','epsc')

figure
% subplot(1,3,2)
title('1D Inviscid Advective N = 41 @ 11xFTT')
plot(Div2.X(Index),Div2.P(Index)/P_0,'LineWidth',2, 'LineStyle','-')
hold on
plot(Div_F4.X(Index),Div_F4.P(Index)/P_0,'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
% plot(DivFilt_G.X(Index),DivFilt_G.P(Index)/rho_0,'LineWidth',2,'LineStyle','-.')
% plot(Div_F2.X(Index),Div_F2.P(Index)/rho_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Hybrid.X(Index),Hybrid.P(Index)/P_0,'LineWidth',2,'LineStyle','-.','color',[0.3010 0.7450 0.9330])
plot(WENO5_X(1:2:end),WENO5_P(1:2:end)/P_0,'--o','LineWidth',2,'MarkerSize',3,'color',[0.6350 0.0780 0.1840])

legend('Div','Div w/F4','Hybrid','WENO5')
legend('Location','best','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{P}/{{P}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
saveas(gca,'Results/1D_Adv_Shima/Ideal/1D_Shima_Ideal_P_2','epsc')

figure
% subplot(1,3,3)
figure
% subplot(1,3,1)
title('1D Inviscid Advective N = 41 @ 11xFTT')
plot(Div2.X(Index),Div2.u(Index)/u_0,'LineWidth',2, 'LineStyle','-')
hold on
plot(Div_F4.X(Index),Div_F4.u(Index)/u_0,'LineWidth',2,'LineStyle',':','color',[0.4660, 0.6740, 0.1880])
% plot(DivFilt_G.X(Index),DivFilt_G.u(Index)/u_0,'LineWidth',2,'LineStyle','-.')
% plot(Div_F2.X(Index),Div_F2.u(Index)/u_0,'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
plot(Hybrid.X(Index),Hybrid.u(Index)/u_0,'LineWidth',2,'LineStyle','-.','color',[0.3010 0.7450 0.9330])
plot(WENO5_X(1:2:end),WENO5_u(1:2:end)/u_0,'--o','LineWidth',2,'MarkerSize',3,'color',[0.6350 0.0780 0.1840])

legend('Div','Div w/F4','Hybrid','WENO5')
legend('Location','best','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{\rho}/{{\rho}_{0}}}$','interpreter','latex')
legend('Location','southwest','box','off')
grid on
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{u}/{{u}_{0}}}$','interpreter','latex')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
saveas(gca,'Results/1D_Adv_Shima/Ideal/1D_Shima_Ideal_u_2','epsc')
