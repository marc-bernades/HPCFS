%% Postprocess 2D Mixing Layer
% clear all; clc; close all;
%% 1D Advective test
% Load results
DATA             = readtable('output_data_2D_MixingLayer_256x128x1_KGP_P_0.066667.csv');
% DATA             = readtable('output_data_2D_MixingLayer_Periodic_Real_20x20x1_KGP_P_3.4558.csv');
% DATA             = readtable('output_data_2D_MixingLayer_64x32x1_KGP_Pt_0.0041667.csv');
DATA             = readtable('output_data_2D_MixingLayer_64x32x1_KGP_et_CR_0.0065148.csv');
DATA             = readtable('output_data_2D_MixingLayer_64x32x1_KGP_Et_CR_0.005718.csv');

% Plot features
mesh      = '256x128';
mesh      = '64x32';
Substance = 'N2';
lambda    = 1/3; %2pi             % Wave length input perturbation
delta_U   = 10; %2;              % Delta speed input field
t_c       = lambda/delta_U;  % Normalized eddy time around
t_end     = 0.1715;   %1.1;               % Normalized cycle time to plot
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Substance);
c_p_Pb    =  3.8591e+03; % Cp at 2xPc (preimposed Pb)
% Index unique Z (2D - X,Y)
Index_Z = unique(DATA.Z); Index_Z = find(DATA.Z == Index_Z(2));
Index = Index_Z;
x = unique(DATA.X);
y = unique(DATA.Y);
[X,Y] = meshgrid(x,y);

% Rho map
z = DATA.rho(Index)/rho_c;
Z = reshape(z,length(y),length(x));
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
colorbar;
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','\rho','/','\rho{}_{c}$'],'Interpreter','latex','fontsize',10);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{Y}/{{L}}}$','interpreter','latex')
set(gca,'linewidth',1)
% set(gca,'fontsize',12)
% ylim([0 1.2])
pbaspect([1.8 1 1])
set(gca, 'OuterPosition', [0,0,1,1])
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_rho'),'epsc')
exportgraphics(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_rho.png'),'Resolution',300)

% T map
z = DATA.T(Index)/T_c;
Z = reshape(z,length(y),length(x));
f = figure; [c,h]=contourf(X,Y,Z,100);
title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
colorbar('Ticks',[0.75 1 1.25 1.5])
clim([0.75 1.5])
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','T','/','T{}_{c}$'],'Interpreter','latex','fontsize',10);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{Y}/{{L}}}$','interpreter','latex')
set(gca,'linewidth',1)
% set(gca,'fontsize',12)
% ylim([0 1.2])
pbaspect([1.8 1 1])
saveas(gca,'Results/2D_MixingLayer/2D_MixingLayer_T','epsc')

% P map
z = DATA.P(Index)/P_c;
Z = reshape(z,length(y),length(x));
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
colorbar;
% clim([0.75 1.5])
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','P','/','P{}_{c}$'],'Interpreter','latex','fontsize',10);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{Y}/{{L}}}$','interpreter','latex')
set(gca,'linewidth',1)
% set(gca,'fontsize',12)
% ylim([0 1.2])
pbaspect([1.8 1 1])
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_P_KGP-et-CR'),'epsc')

% Velocity - Streamwise
z = DATA.u(Index)/delta_U;
Z = reshape(z,length(y),length(x));
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
colorbar;
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','u','/','u{}_{0}$'],'Interpreter','latex','fontsize',10);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{Y}/{{L}}}$','interpreter','latex')
set(gca,'linewidth',1)
% set(gca,'fontsize',12)
% ylim([0 1.2])
pbaspect([1.8 1 1])
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_u_KGP-et-CR'),'epsc')

% Velocity - vertical
z = DATA.v(Index);
Z = reshape(z,length(y),length(x));
f = figure; [c,h]=contourf(X,Y,Z,100);
title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
colorbar;
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$v$'],'Interpreter','latex','fontsize',10);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{Y}/{{L}}}$','interpreter','latex')
set(gca,'linewidth',1)
% set(gca,'fontsize',12)
% ylim([0 1.2])
pbaspect([1.8 1 1])
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_v'),'epsc')

% cp map
z = DATA.c_p(Index)/c_p_Pb;
Z = reshape(z,length(y),length(x));
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
colorbar;
% clim([0.75 1.5])
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','c_p','/','c{}_{p_{Pb}}$'],'Interpreter','latex','fontsize',10);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{Y}/{{L}}}$','interpreter','latex')
set(gca,'linewidth',1)
% set(gca,'fontsize',12)
% ylim([0 1.2])
pbaspect([1.8 1 1])
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_cp'),'epsc')
exportgraphics(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_cp.png'),'Resolution',300)


% sos map
z = DATA.sos(Index);
Z = reshape(z,length(y),length(x));
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
set(h, 'edgecolor','none'); 
colorbar;
% clim([0.75 1.5])
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','c_p','/','c{}_{p_{Pb}}$'],'Interpreter','latex','fontsize',10);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel('${{Y}/{{L}}}$','interpreter','latex')
set(gca,'linewidth',1)
% set(gca,'fontsize',12)
% ylim([0 1.2])
pbaspect([1.8 1 1])
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_cp'),'epsc')
exportgraphics(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , num2str(t_end), 'tc_cp.png'),'Resolution',300)

%% Invariants convergence comparison
DATA_time_64x32        = readtable('output_data_2D_MixingLayer_64x32x1_KGP_P_0.066667_Time.csv');
DATA_time_128x64       = readtable('output_data_2D_MixingLayer_128x64x1_KGP_P_0.066667_Time.csv');
DATA_time_256x128      = readtable('output_data_2D_MixingLayer_256x128x1_KGP_P_0.066667_Time.csv');
% DATA_time_64x32        = readtable('output_data_2D_MixingLayer_Periodic_Real_20x20x1_KGP_P_3.4558_Time.csv');
c_label = {'64x32' ,'128x64','256x128'};

% Subplot figures
figure()
subplot(2,2,1)
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rhoE_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rhoE_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rhoE_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
% legend(c_label)
pbaspect([2.2 1 1])
% legend('Location','southwest','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{E} \rangle}$','interpreter','latex')
% saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhoE','epsc')

subplot(2,2,2)
% figure
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.P_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.P_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.P_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
% legend(c_label)
pbaspect([2.2 1 1])
% legend('Location','northeast','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {P} \rangle}$','interpreter','latex')
% saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_P','epsc')

subplot(2,2,3)
% figure
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rhosPR_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rhosPR_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rhosPR_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
% legend(c_label)
pbaspect([2.2 1 1])
% legend('Location','northeast','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{s} \rangle}$','interpreter','latex')
% saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhoS','epsc')

subplot(2,2,4)
% figure
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rhoui2_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rhoui2_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rhoui2_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
% legend(c_label)
pbaspect([2.2 1 1])
% legend('Location','northeast','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{{u_i}^2} \rangle}$','interpreter','latex')
% saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhou2','epsc')


%% Plot invariants for JCP
figure()
% subplot(2,2,1)
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rho_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rho_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rho_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
legend(c_label)
% pbaspect([2.2 1 1])
legend('Location','northwest','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex','fontsize',16)
ylabel('$\overline{\langle {\rho} \rangle}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rho','epsc')

figure()
% subplot(2,2,1)
plot(DATA_time_64x32.t_vec/t_c, (DATA_time_64x32.rhou_bar -  DATA_time_64x32.rhou_bar(1))/abs(DATA_time_64x32.rhou_bar(1)),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c, (DATA_time_128x64.rhou_bar -  DATA_time_128x64.rhou_bar(1))/abs(DATA_time_128x64.rhou_bar(1)),'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, (DATA_time_256x128.rhou_bar -  DATA_time_256x128.rhou_bar(1))/abs(DATA_time_256x128.rhou_bar(1)), 'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
legend(c_label)
% pbaspect([2.2 1 1])
legend('Location','southeast','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
ylim([-1E-14 2.01E-15])
xlabel('${{t}/{t_c}}$','interpreter','latex','fontsize',16)
ylabel('$\overline{\langle {\rho}{u} \rangle}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhou','epsc')

figure()
% subplot(2,2,1)
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rhoE_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rhoE_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rhoE_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
legend(c_label)
% pbaspect([2.2 1 1])
legend('Location','southwest','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex','fontsize',16)
ylabel('$\overline{\langle {\rho}{E} \rangle}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhoE','epsc')

figure()
% subplot(2,2,1)
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rhoe_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rhoe_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rhoe_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
legend(c_label)
% pbaspect([2.2 1 1])
legend('Location','southwest','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex','fontsize',16)
ylabel('$\overline{\langle {\rho}{e} \rangle}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhoet','epsc')

figure()
% subplot(2,2,1)
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rhoui2_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rhoui2_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rhoui2_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
legend(c_label)
% pbaspect([2.2 1 1])
legend('Location','northeast','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex','fontsize',16)
ylabel('$\overline{\langle {\rho}{{u_i}^2} \rangle}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhou2','epsc')


figure()
% subplot(2,2,1)
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rhos_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rhos_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rhos_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
legend(c_label)
% pbaspect([2.2 1 1])
legend('Location','southwest','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex','fontsize',16)
ylabel('$\overline{\langle {\rho}{s} \rangle}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhos','epsc')

figure()
% subplot(2,2,1)
plot(DATA_time_64x32.t_vec/t_c, DATA_time_64x32.rhosPR_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(DATA_time_128x64.t_vec/t_c,  DATA_time_128x64.rhosPR_norm,'LineWidth',2,'LineStyle','-.','color',[0.4660, 0.6740, 0.1880])
plot(DATA_time_256x128.t_vec/t_c, DATA_time_256x128.rhosPR_norm,'LineWidth',2,'LineStyle',':','color',[0.6350 0.0780 0.1840])
xlim([0 round(t_end)])
legend(c_label)
% pbaspect([2.2 1 1])
legend('Location','northeast','box','off')
set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{t}/{t_c}}$','interpreter','latex','fontsize',16)
ylabel('$\overline{\langle {\rho}{s} \rangle}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\2D_MixingLayer\2D_MixingLayer_Invariants_Conv_rhosPR','epsc')



