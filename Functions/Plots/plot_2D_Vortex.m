function plot_2D_Vortex(u,v,w,X,Y,Z,L_x,L_y,L_z,P, num_grid_x, num_grid_y, num_grid_z, u_exact, v_exact, Fluid)

% Total velocity
u_field = sqrt(u.^2 + v.^2 + w.^2);


% Inner points
u_field = u_field(2:end-1,2:end-1,2:end-1)/Fluid.U_0;
u_plot  = u(2:end-1,2:end-1,2:end-1)/Fluid.U_0;
v_plot  = v(2:end-1,2:end-1,2:end-1)/Fluid.U_0;
x_plot  = X(2:end-1,2:end-1,2:end-1)/Fluid.r_v;
y_plot  = Y(2:end-1,2:end-1,2:end-1)/Fluid.r_v;
p_plot  = P(2:end-1,2:end-1,2:end-1)/Fluid.P_0;
figure()
contourf(x_plot, y_plot,p_plot,5), hold on;
quiver(x_plot,y_plot,u_plot,v_plot)
colorbar
xlabel('x/r_v')
ylabel('y/r_v')
title("2D Vortex Normalized Pressure field (XY) " + num2str(num_grid_x) + "x" + num2str(num_grid_y) + "x" +num2str(num_grid_z))

figure()
contourf(x_plot, y_plot,u_field,5), hold on;
quiver(x_plot,y_plot,u_plot,v_plot)
colorbar
xlabel('x/r_v')
ylabel('y/r_v')
title("2D Vortex Normalized Velocity field (XY) " + num2str(num_grid_x) + "x" + num2str(num_grid_y) + "x" +num2str(num_grid_z))


% Exact solution
% U velocity
figure
x_position_plot = round(num_grid_x/2)+1;
y_position_plot = round(num_grid_y/2)+1;
subplot(2,2,1)
plot(Y(2:end-1,x_position_plot,2:end-1)/Fluid.r_v,u(2:end-1,x_position_plot,2:end-1)/Fluid.U_0,'o'); hold on
plot(Y(2:end-1,x_position_plot,2:end-1)/Fluid.r_v,u_exact(2:end-1,x_position_plot,2:end-1)/Fluid.U_0);
legend(["2D Vortex " + num2str(num_grid_x) + "x" + num2str(num_grid_y)],'Exact')
xlabel('y/r_v')
ylabel('u/U_0 at x/Lx ~ 0.5')
grid on
% xlim([0 1])  
subplot(2,2,2)
plot(X(y_position_plot,2:end-1,2:end-1)/Fluid.r_v,u(y_position_plot,2:end-1,2:end-1)/Fluid.U_0,'o'); hold on
plot(X(y_position_plot,2:end-1,2:end-1)/Fluid.r_v,u_exact(y_position_plot,2:end-1,2:end-1)/Fluid.U_0);
legend(["2D Vortex " + num2str(num_grid_x) + "x" + num2str(num_grid_y)],'Exact')
xlabel('x/r_v')
ylabel('u/u_0 at y/Ly ~ 0.5')
grid on
% xlim([0 1])  

% V velocity
subplot(2,2,3)
plot(Y(2:end-1,x_position_plot,2:end-1)/Fluid.r_v,v(2:end-1,x_position_plot,2:end-1)/Fluid.U_0,'o'); hold on
plot(Y(2:end-1,x_position_plot,2:end-1)/Fluid.r_v,v_exact(2:end-1,x_position_plot,2:end-1)/Fluid.U_0);
legend(["2D Vortex " + num2str(num_grid_x) + "x" + num2str(num_grid_y)],'Exact')
xlabel('y/r_v')
ylabel('v/U_0 at x/Lx ~ 0.5')
grid on
% xlim([0 1])  
subplot(2,2,4)
plot(X(y_position_plot,2:end-1,2:end-1)/Fluid.r_v,v(y_position_plot,2:end-1,2:end-1)/Fluid.U_0,'o'); hold on
plot(X(y_position_plot,2:end-1,2:end-1)/Fluid.r_v,v_exact(y_position_plot,2:end-1,2:end-1)/Fluid.U_0);
legend(["2D Vortex " + num2str(num_grid_x) + "x" + num2str(num_grid_y)],'Exact')
xlabel('x/r_v')
ylabel('v/u_0 at y/Ly ~ 0.5')
grid on
% xlim([0 1])  


end