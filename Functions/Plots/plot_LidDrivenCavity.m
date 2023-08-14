function plot_LidDrivenCavity(u,v,w,X,Y,Z,L_x,L_y,L_z,U_0, num_grid_x, num_grid_y, num_grid_z)

% Inner points normalized
u_plot  = u(2:end-1,2:end-1,2:end-1)/U_0;
v_plot  = v(2:end-1,2:end-1,2:end-1)/U_0;
x_plot  = X(2:end-1,2:end-1,2:end-1)/L_x;
y_plot  = Y(2:end-1,2:end-1,2:end-1)/L_y;

% Total velocity 2D
u_field = sqrt(u_plot.^2 + v_plot.^2);

%% 3D Plot
% Create Countour plot
figure()
contourf(x_plot, y_plot,u_field,20), hold on;
quiver(x_plot,y_plot,u_plot,v_plot)
colorbar

xlabel('x/Lx')
ylabel('y/Ly')
title("Lid-Driven Cavity Normalized Velocity field (XY) " + num2str(num_grid_x) + "x" + num2str(num_grid_y))

%% Reference solution
figure()
title("Ghia et al. reference validation")

subplot(2,1,1)
% Load u vs y results at x/L = 0.5
u_GHIA = csvread("Results/2D_LidDrivenCavity/ghia_u_velocity_solution.csv",2);  
plot(u_GHIA(:,1),u_GHIA(:,2),'ro'); hold on
% Find idx when x/L = 0.5
[val,idx]=min(abs(x_plot(1,:)-0.5));
plot(y_plot(:,idx),u_plot(:,idx),'b')
% [val,idx]=min(abs(y_plot(:,1)-0.5));
% plot(x_plot(idx,:),flip(v_plot(idx,:)),'b')
xlabel('y/L')
ylabel('u/U_0 at x/L = 0.5')
legend('Ghia et al. Re = 100', '3D Compressible Solver')

% Find y = 0.5
v_GHIA = csvread("Results/2D_LidDrivenCavity/ghia_v_velocity_solution.csv",2); 
subplot(2,1,2)
plot(v_GHIA(:,1),v_GHIA(:,2),'ro'); hold on
% Find idx when y/L = 0.5
[val,idx]=min(abs(y_plot(:,1)-0.5));
plot(x_plot(idx,:),v_plot(idx,:),'b')
% [val,idx]=min(abs(x_plot(1,:)-0.5));
% plot(y_plot(:,idx),-u_plot(:,idx),'b')
xlabel('x/L')
ylabel('v/U_0 at y/L = 0.5')
legend('Ghia et al. Re = 100', '3D Compressible Solver')

end