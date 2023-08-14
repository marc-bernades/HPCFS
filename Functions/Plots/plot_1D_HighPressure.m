function plot_1D_HighPressure(rho,c_p,T,e,sos,mu,kappa)

% Load reference data
NIST = csvread("Analytical_solution/nist_nitrogen_4MPa_0-1000K.csv",2);

% Figure comparing main parameters
figure
subplot(2,3,1)
plot(NIST(:,1),NIST(:,3),'o')
hold on; grid on;
plot(T(2,:,2),rho(2,:,2))
xlim([100 200])
ylabel('\rho [kg/m^3]')
xlabel('T [K]')
legend('NIST','Model')

subplot(2,3,2)
plot(NIST(:,1),NIST(:,5),'o')
hold on; grid on;
plot(T(2,:,2),e(2,:,2)/1000 + (387.676 - 76.80))
xlim([100 200])
ylabel('e [kJ/kg]')
xlabel('T [K]')
legend('NIST','Model')

subplot(2,3,3)
plot(NIST(:,1),NIST(:,10),'o')
hold on; grid on;
plot(T(2,:,2),sos(2,:,2))
xlim([100 200])
ylabel('sos [m/s]')
xlabel('T [K]')
legend('NIST','Model')

subplot(2,3,4)
plot(NIST(:,1),NIST(:,9)*1000,'o')
hold on; grid on;
plot(T(2,:,2),c_p(2,:,2))
xlim([100 200])
ylabel('c_p [J/(kg·K)]')
xlabel('T [K]')
legend('NIST','Model')

subplot(2,3,5)
plot(NIST(:,1),NIST(:,12),'o')
hold on; grid on;
plot(T(2,:,2),mu(2,:,2))
xlim([100 200])
ylabel('\mu [Pa·s]')
xlabel('T [K]')
legend('NIST','Model')

subplot(2,3,6)
plot(NIST(:,1),NIST(:,13),'o')
hold on; grid on;
plot(T(2,:,2),kappa(2,:,2))
xlim([100 200])
ylabel('\kappa [W/(m·K)]')
xlabel('T [K]')
legend('NIST','Model')

end