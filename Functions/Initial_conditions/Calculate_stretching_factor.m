function A = Calculate_stretching_factor(y_plus, u_tau, Re_tau, delta,num_grid, X_0, L)

% Calculate parameters
nu = u_tau*delta/Re_tau;
y = y_plus*nu/u_tau;

l = 1; %First grid point
eta = (l-0.5)/num_grid;
A = (y - X_0 - 2*delta*eta)./((0.5*L - L*eta)*(1-eta)*eta);

end