function plot_3D_TGV_Invariants_KGP(t_vec, rho_norm, rhou_bar, rhoE_norm, P_norm, rhos_norm, rhoui2_norm, ...
    rhoE2_norm, rhoe2_norm, rhos2_norm, varargin)

% Legend label
c_label = varargin{end-1};

% xLimit t/t_C
t_c   = varargin{end-2};
t_max = varargin{end};

%% FULL INVARIANT PLOTS
figure()
% t_c = 2*pi;
subplot(3,3,1)
plot(t_vec/t_c, rho_norm)
hold on
plot(varargin{1}/t_c,varargin{2})
%         xlim([0 max(t_vec)])
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho} \rangle}$','interpreter','latex')

subplot(3,3,2)
plot(t_vec/t_c, rhou_bar)
hold on
plot(varargin{1}/t_c,varargin{3})
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{{\rho}{u}}$','interpreter','latex')

subplot(3,3,3)
plot(t_vec/t_c, rhoE_norm)
hold on
plot(varargin{1}/t_c,varargin{4})
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{E} \rangle}$','interpreter','latex')

subplot(3,3,4)
plot(t_vec/t_c, P_norm)
hold on
plot(varargin{1}/t_c,varargin{5})
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {P} \rangle}$','interpreter','latex')

subplot(3,3,5)
plot(t_vec/t_c, rhos_norm)
hold on
plot(varargin{1}/t_c,varargin{6})
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{s} \rangle}$','interpreter','latex')

subplot(3,3,6)
plot(t_vec/t_c, rhoui2_norm)
hold on
plot(varargin{1}/t_c,varargin{7})
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{{u_i}^2} \rangle}$','interpreter','latex')

subplot(3,3,7)
plot(t_vec/t_c, rhoE2_norm)
hold on
plot(varargin{1}/t_c,varargin{8})
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{{E}^2} \rangle}$','interpreter','latex')

subplot(3,3,8)
plot(t_vec/t_c, rhoe2_norm)
hold on
plot(varargin{1}/t_c,varargin{9})
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{{e}^2} \rangle}$','interpreter','latex')

subplot(3,3,9)
plot(t_vec/t_c, rhos2_norm)
hold on
plot(varargin{1}/t_c,varargin{10})
xlim([0 t_max])
legend(c_label)
xlabel('${{t}/{t_c}}$','interpreter','latex')
ylabel('$\overline{\langle {\rho}{{s}^2} \rangle}$','interpreter','latex')

%         pbaspect([1.8 1 1])
% saveas(gca,'Results\3D_TGV\t_20FTT\3D_TGV_Invariants_KGP','epsc')

%% Reduced invariants plot
figure()
% t_c = 2*pi;
% subplot(2,2,1)
plot(t_vec/t_c, rhoE_norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(varargin{1}/t_c,varargin{4},'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
xlim([0 t_max])
legend(c_label)
legend('Location','northwest','box','off')
pbaspect([1.8 1 1])
xlabel('${{t}/{t_c}}$','interpreter','latex','FontWeight','bold')
ylabel('$\overline{\langle {\rho}{E} \rangle}$','interpreter','latex','FontWeight','bold')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,'Results\3D_TGV\t_20FTT\3D_TGV_Invariants_KGP_reduced_E','epsc')

% subplot(2,2,2)
figure
plot(t_vec/t_c, P_norm, 'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(varargin{1}/t_c,varargin{5},'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
xlim([0 t_max])
legend(c_label)
legend('Location','northwest','box','off')
pbaspect([1.8 1 1])
xlabel('${{t}/{t_c}}$','interpreter','latex','FontWeight','bold')
ylabel('$\overline{\langle {P} \rangle}$','interpreter','latex','FontWeight','bold')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,'Results\3D_TGV\t_20FTT\3D_TGV_Invariants_KGP_reduced_P','epsc')

% subplot(2,2,3)
figure
plot(t_vec/t_c, rhos_norm, 'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(varargin{1}/t_c,varargin{6},'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
xlim([0 t_max])
legend(c_label)
legend('Location','southwest','box','off')
pbaspect([1.8 1 1])
xlabel('${{t}/{t_c}}$','interpreter','latex','FontWeight','bold')
ylabel('$\overline{\langle {\rho}{s} \rangle}$','interpreter','latex','FontWeight','bold')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,'Results\3D_TGV\t_20FTT\3D_TGV_Invariants_KGP_reduced_s','epsc')

% subplot(2,2,4)
figure
plot(t_vec/t_c, rhoui2_norm, 'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
hold on
plot(varargin{1}/t_c,varargin{7},'LineWidth',2,'LineStyle','--','color',[0.3010 0.7450 0.9330])
xlim([0 t_max])
legend(c_label)
legend('Location','southwest','box','off')
pbaspect([1.8 1 1])
xlabel('${{t}/{t_c}}$','interpreter','latex','FontWeight','bold')
ylabel('$\overline{\langle {\rho}{{u_i}^2} \rangle}$','interpreter','latex','FontWeight','bold')
set(gca,'linewidth',1.5)
set(gca,'fontsize',14)
saveas(gca,'Results\3D_TGV\t_20FTT\3D_TGV_Invariants_KGP_reduced_ke','epsc')


%         pbaspect([1.8 1 1])
% saveas(gca,'Results\3D_TGV\t_20FTT\3D_TGV_Invariants_KGP_reduced','epsc')


end