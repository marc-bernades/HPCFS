%% Postprocess 2D Mixing Layer
% clear all; clc; close all;
%% 1D Advective test
% Load results
DATA1             = readtable('output_data_2D_MixingLayer_64x32x1_KGP_Pt_0.0066667.csv');
DATA2             = readtable('output_data_2D_MixingLayer_64x32x1_KGP_Et_0.0066667.csv');

% Plot features
DATA      = DATA1;
scheme    = 'KGP-Pt';
mesh      = '64x32'; %'256x128';
Substance = 'N2';
lambda    = 1/3; %2pi             % Wave length input perturbation
delta_U   = 10; %2;              % Delta speed input field
t_c       = lambda/delta_U;  % Normalized eddy time around
t_end     = 0.004166; %2   1.1;               % Normalized cycle time to plot
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
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , scheme, '_', num2str(t_end), 'tc_rho'),'epsc')
exportgraphics(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , scheme, '_',  num2str(t_end), 'tc_rho.png'),'Resolution',300)

% T map
z = DATA.T(Index)/T_c;
Z = reshape(z,length(y),length(x));
f = figure; [c,h]=contourf(X,Y,Z,100);
% title('${{N_2},~{P}/{P_c}={2},~{T}/{T_c}\equiv{0.75-1.5},~{256}x{128},~{CFL=0.3},~{t/{t_c}}{=}~{2}}$','Interpreter','latex','fontsize',10)
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
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , scheme, '_', num2str(t_end), 'tc_T'),'epsc')
exportgraphics(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , scheme, '_',  num2str(t_end), 'tc_T.png'),'Resolution',300)

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
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , scheme, '_', num2str(t_end), 'tc_P'),'epsc')
exportgraphics(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , scheme, '_',  num2str(t_end), 'tc_P.png'),'Resolution',300)

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
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , scheme, '_', num2str(t_end), 'tc_u'),'epsc')
exportgraphics(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , scheme, '_',  num2str(t_end), 'tc_u.png'),'Resolution',300)


%% Pressure comparisons
clabel    = {'KGP-Pt','KGP-Et'}; % DATA 1, DATA 2

% Plots at X/L = 0
target = 0;

% Normalized pressure
figure
z1 = DATA1.P(Index)/P_c;
Z1 = reshape(z1,length(y),length(x));
z2 = DATA2.P(Index)/P_c;
Z2 = reshape(z2,length(y),length(x));
[ d, ix ] = min( abs( x-target ) );
plot(Y(:,ix),Z1(:,ix),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410]); hold on
plot(Y(:,ix),Z2(:,ix),'LineWidth',2,'LineStyle','-','color',[0.6350 0.0780 0.1840]);
xlabel('${{Y}/{L}}$','interpreter','latex')
ylabel(['$','P','/','P{}_{c}$'],'interpreter','latex')
xlim([-0.25 0.25])
pbaspect([1.8 1 1])
legend(clabel)
legend('Location','best','box','off')
% grid on
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , 'P_y'),'epsc')

% Normalized density
figure
z1 = DATA1.rho(Index)/rho_c;
Z1 = reshape(z1,length(y),length(x));
z2 = DATA2.rho(Index)/rho_c;
Z2 = reshape(z2,length(y),length(x));
[ d, ix ] = min( abs( x-target ) );
plot(Y(:,ix),Z1(:,ix),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410]); hold on
plot(Y(:,ix),Z2(:,ix),'LineWidth',2,'LineStyle','-','color',[0.6350 0.0780 0.1840]);
xlabel('${{Y}/{L}}$','interpreter','latex')
ylabel(['$','\rho','/','\rho{}_{c}$'],'interpreter','latex')
xlim([-0.25 0.25])
pbaspect([1.8 1 1])
legend(clabel)
legend('Location','best','box','off')
% grid on
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

axes('position',[.59 .3 .28 .25])
plot(Y(:,ix),Z1(:,ix),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410]); hold on
plot(Y(:,ix),Z2(:,ix),'LineWidth',2,'LineStyle','-','color',[0.6350 0.0780 0.1840]);
xlim([0.15, 0.25])
ylim([2.66,2.7])
% pbaspect([1.8 1 1])
% set(gca,'linewidth',1)
% set(gca,'fontsize',10)
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , 'rho_y'),'epsc')



% Plots at Y/L = 0.125
target = 0.125;
figure
z1 = DATA1.P(Index)/P_c;
Z1 = reshape(z1,length(y),length(x));
z2 = DATA2.P(Index)/P_c;
Z2 = reshape(z2,length(y),length(x));
[ d, ix ] = min( abs( y-target ) );
plot(X(ix,:),Z1(ix,:),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410]); hold on
plot(X(ix,:),Z2(ix,:),'LineWidth',2,'LineStyle','-','color',[0.6350 0.0780 0.1840]);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel(['$','P','/','P{}_{c}$'],'interpreter','latex')
xlim([-0.5 0.5])
pbaspect([1.8 1 1])
legend(clabel)
legend('Location','best','box','off')
% grid on
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , 'P_x'),'epsc')

figure
z1 = DATA1.rho(Index)/rho_c;
Z1 = reshape(z1,length(y),length(x));
z2 = DATA2.rho(Index)/rho_c;
Z2 = reshape(z2,length(y),length(x));
[ d, ix ] = min( abs( y-target ) );
plot(X(ix,:),Z1(ix,:),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410]); hold on
plot(X(ix,:),Z2(ix,:),'LineWidth',2,'LineStyle','-','color',[0.6350 0.0780 0.1840]);
xlabel('${{X}/{L}}$','interpreter','latex')
ylabel(['$','\rho','/','\rho{}_{c}$'],'interpreter','latex')
xlim([-0.5 0.5])
pbaspect([1.8 1 1])
legend(clabel)
legend('Location','best','box','off')
% grid on
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,strcat('Results/2D_MixingLayer/2D_MixingLayer_',mesh, '_' , 'rho_x'),'epsc')







