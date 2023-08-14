figure()
ny = 10;
title('Evaluation enthalpy splitting - 2D TGV at 0.1s')
subplot(2,2,1)
title('Evaluation enthalpy splitting - 2D TGV at 0.1s')

plot(X(ny,:,2)/L_x, u_ref_KGP(ny,:,2)/Fluid.U_0,'k','LineWidth',5); hold on
plot(X(ny,:,2)/L_x,u_ref_F(ny,:,2)/Fluid.U_0,'r','LineWidth',3)
plot(X(ny,:,2)/L_x,u_ref_C(ny,:,2)/Fluid.U_0,'b','LineWidth',1.5)
ylabel('u/U_0')
xlabel('x/L_x')
legend('KGP','F','C')

subplot(2,2,2)
plot(X(ny,:,2)/L_x, E_ref_KGP(ny,:,2),'k','LineWidth',5); hold on
plot(X(ny,:,2)/L_x,E_ref_F(ny,:,2),'r','LineWidth',3)
plot(X(ny,:,2)/L_x,E_ref_C(ny,:,2),'b','LineWidth',1.5)
ylabel('E [Pa]')
xlabel('x/L_x')
legend('KGP','F','C')

subplot(2,2,3)
plot(X(ny,:,2)/L_x, u_ref_KGP(ny,:,2)/Fluid.U_0,'k','LineWidth',3); hold on
plot(X(ny,:,2)/L_x,u_enth_KGP(ny,:,2)/Fluid.U_0,'r','LineWidth',1.5)
ylabel('u/U_0')
xlabel('x/L_x')
legend('E KGP','h split KGP')


subplot(2,2,4)
plot(X(ny,:,2)/L_x, E_ref_KGP(ny,:,2),'k','LineWidth',3); hold on
plot(X(ny,:,2)/L_x,E_enth_KGP(ny,:,2),'r','LineWidth',1.5)
ylabel('E [Pa]')
xlabel('x/L_x')
legend('E KGP','h split KGP')


