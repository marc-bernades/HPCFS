function plot_3D_TGV(u,v,w,X,Y,Z,L_x,L_y,L_z,U_0, num_grid_x, num_grid_y, num_grid_z, z_plane, Fluid, P)

% Total velocity
u_field = sqrt(u.^2 + v.^2 + w.^2);

% Inner points
u_field = u_field(2:end-1,2:end-1,z_plane);
u_plot  = u(2:end-1,2:end-1,z_plane)/U_0;
v_plot  = v(2:end-1,2:end-1,z_plane)/U_0;
x_plot  = X(2:end-1,2:end-1,z_plane)/L_x;
y_plot  = Y(2:end-1,2:end-1,z_plane)/L_y;
p_plot  = P(2:end-1,2:end-1,z_plane)/Fluid.P_0;

figure()
contourf(x_plot, y_plot,u_field,5), hold on;
quiver(x_plot,y_plot,u_plot,v_plot)
colorbar
xlabel('x/Lx')
ylabel('y/Ly')
title("3D TGV Normalized Velocity field (XY) " + num2str(num_grid_x) + "x" + num2str(num_grid_y) + "x" +num2str(num_grid_z) + " at z plane " + num2str(Z(1,1,z_plane)/pi) + "x pi")

figure()
contourf(x_plot, y_plot,p_plot,5), hold on;
quiver(x_plot,y_plot,u_plot,v_plot)
colorbar
xlabel('x/Lx')
ylabel('y/Ly')
title("3D TGV Normalized Pressure (XY) " + num2str(num_grid_x) + "x" + num2str(num_grid_y) + "x" +num2str(num_grid_z) + " at z plane " + num2str(Z(1,1,z_plane)/pi) + "x pi")




end