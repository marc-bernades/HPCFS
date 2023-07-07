function [grid, rho, u, v, w, E, rhou, rhov, rhow, rhoE, rhoe, P, T, mu, kappa, sos, rhou_0, rhov_0, rhow_0, rhoE_0, rhoe_0, ...
    rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv, rhoe_conv, dP_rhou, dP_rhov, dP_rhow, dP_rhoE, ...
    rhou_vis, rhov_vis, rhow_vis, rhoE_vis, f_rhou, f_rhov, f_rhow, f_rhoE] = Initialize_Variables(num_grid_x, num_grid_y, num_grid_z)

% Meshgrid sets x to work horizontal and y vertical
num_grid_x_0 = num_grid_x;
num_grid_x   = num_grid_y;
num_grid_y   = num_grid_x_0;

% Initialize grid
grid = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

% Initialize primative variables
rho = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
u   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
v   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
w   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
E   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

% Conserved variables
rhou   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhov   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhow   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhoE   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhoe   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

% Thermodynamic variable
P     = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
T     = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
mu    = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
kappa = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
sos   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

% Time integration (current state) variables
rho_0  = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhou_0 = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhov_0 = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhow_0 = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhoE_0 = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhoe_0 = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

% Convective terms
rho_conv  = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhou_conv = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhov_conv = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhow_conv = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhoE_conv = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhoe_conv = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

% Pressure gradient terms
dP_rhou = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
dP_rhov = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
dP_rhow = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
dP_rhoE = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

% Viscous terms
rhou_vis = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhov_vis = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhow_vis = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
rhoE_vis = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

% Source terms
f_rhou   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
f_rhov   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
f_rhow   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);
f_rhoE   = zeros(num_grid_x + 2, num_grid_y + 2, num_grid_z + 2);

end

