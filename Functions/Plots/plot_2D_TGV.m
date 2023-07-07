function plot_2D_TGV(u,v,w,X,Y,Z,L_x,L_y,L_z,U_0, num_grid_x, num_grid_y, num_grid_z, u_exact, v_exact, Fluid,c_p,P,rho,T)

% Total velocity
u_field = sqrt(u.^2 + v.^2 + w.^2);
% [rho,e,ke,E,sos,c_v,c_p, mu, kappa] = Initialize_Thermodynamics(P,T,Fluid.R_specific,Fluid.gamma,u,v,w, mu, Fluid.mu_0, kappa, Fluid.kappa_0);
% Z = P./(rho.*Fluid.R_specific.*T);
% Equi = (c_p - Z.*Fluid.R_specific)./(Z.*Fluid.R_specific);
% Equi = Equi(2:end-1,2:end-1,2:end-1);
% Inner points
u_field = u_field(2:end-1,2:end-1,2:end-1);
u_plot  = u(2:end-1,2:end-1,2:end-1)/U_0;
v_plot  = v(2:end-1,2:end-1,2:end-1)/U_0;
x_plot  = X(2:end-1,2:end-1,2:end-1)/L_x;
y_plot  = Y(2:end-1,2:end-1,2:end-1)/L_y;
figure()
% contourf(x_plot, y_plot,Equi,5), hold on;
contourf(x_plot, y_plot,u_field,5), hold on;
quiver(x_plot,y_plot,u_plot,v_plot)
colorbar
xlabel('x/Lx')
ylabel('y/Ly')
title("3D TGV Normalized Velocity field (XY) " + num2str(num_grid_x) + "x" + num2str(num_grid_y) + "x" +num2str(num_grid_z))


% Exact solution
% U velocity
figure
x_position_plot = round(num_grid_x/2)+1;
y_position_plot = round(num_grid_y/2)+1;
subplot(2,2,1)
plot(Y(2:end-1,x_position_plot,2:end-1)/L_y,u(2:end-1,x_position_plot,2:end-1)/Fluid.U_0,'o'); hold on
plot(Y(2:end-1,x_position_plot,2:end-1)/L_y,u_exact(2:end-1,x_position_plot,2:end-1)/Fluid.U_0);
legend(["2D TGV " + num2str(num_grid_x) + "x" + num2str(num_grid_y)],'Exact')
xlabel('y/L_y')
ylabel('u/U_0 at x/Lx ~ 0.5')
grid on
xlim([0 1])  
subplot(2,2,2)
plot(X(y_position_plot,2:end-1,2:end-1)/L_x,u(y_position_plot,2:end-1,2:end-1)/Fluid.U_0,'o'); hold on
plot(X(y_position_plot,2:end-1,2:end-1)/L_x,u_exact(y_position_plot,2:end-1,2:end-1)/Fluid.U_0);
legend(["2D TGV " + num2str(num_grid_x) + "x" + num2str(num_grid_y)],'Exact')
xlabel('x/L_x')
ylabel('u/u_0 at y/Ly ~ 0.5')
grid on
xlim([0 1])  

% V velocity
subplot(2,2,3)
plot(Y(2:end-1,x_position_plot,2:end-1)/L_y,v(2:end-1,x_position_plot,2:end-1)/Fluid.U_0,'o'); hold on
plot(Y(2:end-1,x_position_plot,2:end-1)/L_y,v_exact(2:end-1,x_position_plot,2:end-1)/Fluid.U_0);
legend(["2D TGV " + num2str(num_grid_x) + "x" + num2str(num_grid_y)],'Exact')
xlabel('y/L_y')
ylabel('v/U_0 at x/Lx ~ 0.5')
grid on
xlim([0 1])  
subplot(2,2,4)
plot(X(y_position_plot,2:end-1,2:end-1)/L_x,v(y_position_plot,2:end-1,2:end-1)/Fluid.U_0,'o'); hold on
plot(X(y_position_plot,2:end-1,2:end-1)/L_x,v_exact(y_position_plot,2:end-1,2:end-1)/Fluid.U_0);
legend(["2D TGV " + num2str(num_grid_x) + "x" + num2str(num_grid_y)],'Exact')
xlabel('x/L_x')
ylabel('v/u_0 at y/Ly ~ 0.5')
grid on
xlim([0 1])  


end