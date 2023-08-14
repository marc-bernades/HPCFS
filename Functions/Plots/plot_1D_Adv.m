function plot_1D_Adv(rho,rho_0,X,L_x,U_0, num_grid_x,time,Fluid,Substance)


% Inner points
rho_plot  = rho(2:end-1,2:end-1,2:end-1)/rho_0;
x_plot    = X(2:end-1,2:end-1,2:end-1)/L_x;
[rho_exact, ~, ~, ~, ~, ~, ~, ~, ~ ] = Reference_Solution_1D_Adv(X,Fluid,Substance);
rho_exact = rho_exact(2:end-1,2:end-1,2:end-1)/rho_0;

% Plot
figure()
plot(x_plot, rho_plot,'o'), hold on;
plot(x_plot, rho_exact)
xlabel('x/Lx')
ylabel('rho/rho_0')
title("1D Advective " + num2str(num_grid_x) + "x at time " + num2str(round(time)) + " s total turn over time " + num2str(round(time/(L_x/U_0)),1))
legend('1D Adv','Exact')
end