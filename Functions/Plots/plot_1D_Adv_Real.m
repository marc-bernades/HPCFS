function plot_1D_Adv(rho,u,E,P,T,e,ke,X,L,Substance,Fluid)    
% Update X/L
X = X./L;

% Load reference solution
[rho_ref, P_ref, T_ref, e_ref, u_ref,v_ref,w_ref, ke_ref, E_ref ] = Reference_Solution_1D_Adv(X,Fluid,Substance);

% Overlay validation
figure
subplot(2,3,1)
plot(X(2,:,2),rho_ref(2,:,2),'k-',LineWidth=5)
hold on;
plot(X(2,:,2),rho(2,:,2),LineWidth=2.5)
legend('Ma et. al 2017 Smooth','1D Adv')
xlabel('X/L_x')
ylabel('\rho [kg/m^3]')
grid on;
xlim([0,1])

subplot(2,3,2)
plot(X(2,:,2),u_ref(2,:,2),'k-',LineWidth=5)
hold on;
plot(X(2,:,2),u(2,:,2),LineWidth=2.5)
legend('Ma et. al 2017 Smooth','1D Adv')
xlabel('X/L_x')
ylabel('u [m/s]')
grid on;
xlim([0,1])

subplot(2,3,3)
plot(X(2,:,2),E_ref(2,:,2),'k-',LineWidth=5)
hold on;
plot(X(2,:,2),E(2,:,2),LineWidth=2.5)
legend('Ma et. al 2017 Smooth','1D Adv')
xlabel('X/L_x')
ylabel('E [J]')
grid on;
xlim([0,1])

subplot(2,3,4)
plot(X(2,:,2),P_ref(2,:,2),'k-',LineWidth=5)
hold on;
plot(X(2,:,2),P(2,:,2),LineWidth=2.5)
legend('Ma et. al 2017 Smooth','1D Adv')
xlabel('X/L_x')
ylabel('P [Pa]')
grid on;
xlim([0,1])

subplot(2,3,5)
plot(X(2,:,2),T_ref(2,:,2),'k-',LineWidth=5)
hold on;
plot(X(2,:,2),T(2,:,2),LineWidth=2.5)
legend('Ma et. al 2017 Smooth','1D Adv')
xlabel('X/L_x')
ylabel('T [K]')
grid on;
xlim([0,1])

subplot(2,3,6)
plot(X(2,:,2),e_ref(2,:,2),'k-',LineWidth=5)
hold on;
plot(X(2,:,2),e(2,:,2),LineWidth=2.5)
legend('Ma et. al 2017 Smooth','1D Adv')
xlabel('X/L_x')
ylabel('e [J/kg]')
grid on;
xlim([0,1])

end