%% 3D Compressible solver %%
close all; clear all; clc

%% SET PROBLEM PARAMETERS %%
% Thermodynamics framework
bSolver           = 'Ideal';                    % Define Ideal or Real Gas Model
% Substance
Substance         = 'N2';                      % Fluid selected substance
HP_model          = 'Constant';                % Transport coefficients model: 'Constant', 'LowPressure', 'HighPressure'


%% SET PROBLEM PARAMETERS %%
% Fluid properties
if strcmp(bSolver,'Ideal')
    % IDEAL
    Fluid.R_specific  = 287.058;                % Specific gas constant of fluid [J/(kg K)]
    Fluid.gamma       = 1.4;                    % Heat capacity ratio of the fluid
    Fluid.Re          = 100;                    % Re number bulk
    Fluid.Ma          = 0.01;                    % Mach number
    Fluid.Pr          = 0.71;                   % Pr number
    Fluid.c_p         = Fluid.gamma*Fluid.R_specific/(Fluid.gamma-1);
    Fluid.L           = 1;                      % Reference length [m]
    Fluid.U_0         = 1;                      % Reference velocity [m/s]
    Fluid.rho_0       = 1;                      % Reference density
    Fluid.T_0         = Fluid.U_0/(Fluid.Ma^2*Fluid.gamma*Fluid.R_specific); % Reference temperature [k]
    Fluid.P_0         = Fluid.rho_0*Fluid.R_specific*Fluid.T_0;              % Reference pressure [Pa]
    Fluid.mu_0        = Fluid.rho_0*Fluid.U_0*Fluid.L/Fluid.Re;              % Dynamic viscosity of the fluid [Pa s]
    Fluid.kappa_0     = Fluid.c_p*Fluid.mu_0/Fluid.Pr;                       % Thermal conductivity of the fluid
else
    % REAL

end

% Problem parameters
x_0           = 0;                          % Domain origin in x-direction [m]
y_0           = 0;                          % Domain origin in y-direction [m]
z_0           = 0;                          % Domain origin in z-direction [m]
L_x           = 1.0*Fluid.L;                % Size of domain in x-direction [m]
L_y           = 1.0*Fluid.L;                % Size of domain in y-direction [m]
L_z           = 0.1*Fluid.L;                % Size of domain in z-direction [m]
t_0           = 0.0;                        % Initial time [s]
t_end         = 10;                       % Final time [s]
Test          = 'Lid-Driven_Cavity';        % Test assessment for Initial conditions

% Computational parameters
num_grid_x  = 32;
num_grid_y  = 32;
num_grid_z  = 1;
CFL         = 0.9;                  % CFL time-step stability coefficient
n_iter_max  = 1E10;                 % Maximum number of time iterations
output_iter = 1E2;                  % Output data every output_iter
RK_order    = 4;                    % Time-integration Runge-Kutta order
% Stretching
Fluid.A_x = 0;
Fluid.A_y = 0;
Fluid.A_z = 0;

% Source terms
ST.f_rhou = 0;
ST.f_rhov = 0;
ST.f_rhow = 0;
ST.f_rhoE = 0;

% Define Controller
K.b_active      = 0;    % Boolean activating controller
ST.f_rhou       = 0; % Set initial Source term value

% Pressure spurious oscillation
bPressureModel          = 0;        % Use pressure model to avoid spurious oscillations

% Filter conserved variables
bFilter.Gate            = 0;           % Filter applied (2nd order)  
bFilter.Type            = 'Gaussian';  % Implicit or Gaussian

% Convective scheme: Matrix or flux form
% Matrix form
bEnergySplit            = 1;        % Assess Split Energy Equation
bEnergySplitEnthalpy    = 0;        % Assess Inviscid Energy with Pressure gradient > Enthalpy
% Flux form
bFlux.gate              = 0;        % Inviscid flux wise
bFlux.type              = 'WENO5';    % Type of flux 'HLLC' 'HLLC_plus', 'KGP', 'WENO5'
bFlux.hybrid            = 0;         % Hybrid between energy split and flux

% Convection scheme definition: [alpha (CD), beta (Cfi), gamma(Cu), delta(Crho), epsilon(CL)]
scheme.def      = 1;                   % Scheme: Central = 0, KGP = 1, F = 2, C = 3
scheme.mass     = [1/2 0 1/2 0 0];
scheme.momentum = [1/4 1/4 1/4 1/4 0];
scheme.e        = [1/4 1/4 1/4 1/4 0];
scheme.k        = [1/4 1/4 1/4 1/4 0];
scheme.p        = [1/4 1/4 1/4 1/4 0];

% Viscous Scheme
bViscous_5PointsStencil = 0;        % Gate for 5-points stencil choice
bViscosityVariable      = 1;        % Include effects of variable viscosity over space
bKappaVariable          = 1;        % Include effects of variable thermal conductivity

% Output file name
name_file_out = ['output_data_' Test '_' num2str(num_grid_x) ...
    'x' num2str(num_grid_y) 'x' num2str(num_grid_z) '_KGP_Et'];


% Boundary conditions
BC_types = {'Periodic','Neumann','Dirichlet'}; % Select BC_Types{1},{2},{3}
% Define BC for BC_type & [u,v,w,P,T]
% Dirichlet: Set P<0 for impermeable boundary
% Dirichlet: Set T<=0 for no Temperature effect
BC.west   = {BC_types{3},0,0,0,-1,Fluid.T_0 };
BC.east   = {BC_types{3},0,0,0,-1,Fluid.T_0 };
BC.south  = {BC_types{3},0,0,0,-1,Fluid.T_0 };
BC.north  = {BC_types{3},1,0,0,-1,Fluid.T_0 };
BC.back   = {BC_types{1},0,0,0,0,0};
BC.front  = {BC_types{1},0,0,0,0,0};



%% 3D COMPRESSIBLE SOLVER
[u,v,w,E,rho,mu,kappa,c_v,c_p,P,T,ke,e,sos,Beta_T, Beta_v, Beta_s, Alpha_p,time,t_vec,ke_total,X,Y,Z, Invariants] = Compressible_Solver(Fluid,bSolver,x_0,y_0,z_0,L_x,L_y, L_z,...
        t_0,t_end,name_file_out,num_grid_x,num_grid_y,num_grid_z, CFL, n_iter_max, output_iter, RK_order, scheme, BC, ST, K, ...
        bViscous_5PointsStencil, bViscosityVariable, bKappaVariable, bEnergySplit, bEnergySplitEnthalpy, Test, Substance, HP_model, bPressureModel, bFlux, bFilter);

%% Data output
DataOutput(u,v,w,E,rho,mu,kappa,P,T,ke,e,c_p,c_v,sos,Beta_T,Beta_v,Beta_s,Alpha_p,X,Y,Z,t_vec,ke_total,Invariants,name_file_out);
    
  
% Solution Plots
plot_LidDrivenCavity(u,v,w,X,Y,Z,L_x,L_y,L_z,Fluid.U_0, num_grid_x, num_grid_y, num_grid_z);



