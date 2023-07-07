%% 3D Compressible solver %%
close all; clear all; clc

%% SET PROBLEM PARAMETERS %%
% Thermodynamics framework
bSolver           = 'Real';                    % Define Ideal or Real Gas Model
% Substance
Substance         = 'N2';                      % Fluid selected substance
HP_model          = 'Constant';                % Transport coefficients model: 'Constant', 'LowPressure', 'HighPressure'

% Fluid properties
if strcmp(bSolver,'Ideal')
    % IDEAL
    Fluid.R_specific  = 287.058;                % Specific gas constant of fluid [J/(kg K)]
    Fluid.gamma       = 1.4;                    % Heat capacity ratio of the fluid
    Fluid.Ma          = 0.01;                   % Mach number
    Fluid.U_0         = 1;                      % Reference velocity [m/s]
    Fluid.rho_0       = 2;                      % Reference density
    Fluid.rho_min     = 1;                      % Min density non-linear
    Fluid.rho_max     = 3;                      % Max density non-linear
%     Fluid.T0_min      = 500;                    % Min density non-linear
%     Fluid.T0_max      = 1000;                   % Max density non-linear
    Fluid.P_0         = 1;                    % Reference pressure [Pa]
    Fluid.T_0         = Fluid.P_0/(Fluid.rho_0*Fluid.R_specific);    % Reference temperature [k]
    Fluid.mu_0        = 0.0;                    % Dynamic viscosity of the fluid [Pa s]
    Fluid.kappa_0     = 0.0;                    % Thermal conductivity of the fluid
else
    % REAL
%     Fluid.R_specific  = 287.058; 
%     Fluid.gamma       = 1.4;                  % Heat capacity ratio of the fluid
    Fluid.U_0         = 1;                      % Reference velocity [m/s]
    Fluid.P_0         = 5E6;                    % Rereference Pressure [Pa]
    Fluid.rho_min     = 56.9;                   % Min density non-linear
    Fluid.rho_max     = 793.1;                  % Max density non-linear
%     Fluid.T0_min     = 500;                   % Min density non-linear
%     Fluid.T0_max     = 1000;                  % Max density non-linear
    Fluid.Case        = 'Smooth';               % Non-linear case "Smooth" as per Shima et al.
    Fluid.mu_0        = 0;
    Fluid.kappa_0     = 0;
end
% Problem parameters
x_0           = 0;                          % Domain origin in x-direction [m]
y_0           = 0;                          % Domain origin in y-direction [m]
z_0           = 0;                          % Domain origin in z-direction [m]
L_x           = 1;                          % Size of domain in x-direction [m]
L_y           = 0.01*1;                     % Size of domain in y-direction [m]
L_z           = 0.01*1;                     % Size of domain in z-direction [m]
t_0           = 0.0;                        % Initial time [s]
t_c           = L_x/Fluid.U_0;              % Characteristic time [s]
t_end         = 0.01;                       % Final time [s]
Test          = '1D_Adv';                   % Test assessment for Initial conditions and to only transport rhp

% Computational parameters
num_grid_x  = 41;
num_grid_y  = 1;
num_grid_z  = 1;
CFL         = 0.3;                   % CFL time-step stability coefficient
n_iter_max  = 1E10;                  % Maximum number of time iterations
output_iter = t_c/100:t_c/100:t_end; % Output data every selected time step
RK_order    = 4;                     % Time-integration Runge-Kutta order

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

%% Pressure-based NS - Pressure spurious oscillation
bPressureModel          = 1;        % Use pressure model to avoid spurious oscillations

%% Filter conserved variables
bFilter.Gate            = 0;           % Filter applied (2nd order)  
bFilter.Type            = 'Implicit';  % Implicit or Gaussian

%% Convective scheme
% Convective scheme: Matrix or flux form
% Matrix form
bEnergySplit            = 1;        % Assess Split Energy Equation
bEnergySplitEnthalpy    = 0;        % Assess Inviscid Energy with Pressure gradient > Enthalpy
% Flux form
bFlux.gate              = 0;        % Inviscid flux wise
bFlux.type              = 'HLLC';   % Type of flux 'HLLC_plus', 'KGP', 'WENO5'
bFlux.hybrid            = 0;        % Hybrid between energy split and flux

% Convection scheme definition: [alpha (CD), beta (Cfi), gamma(Cu), delta(Crho), epsilon(CL)]
scheme.def      = 1;                   
scheme.mass     = [1/2 0 1/2 0 0];
scheme.momentum = [1/4 1/4 1/4 1/4 0];
scheme.e        = [1/4 1/4 1/4 1/4 0]; 
scheme.k        = [1/4 1/4 1/4 1/4 0];
scheme.p        = [1/4 1/4 1/4 1/4 0]; 

%% Viscous Scheme
bViscous_5PointsStencil = 0;        % Gate for 5-points stencil choice
bViscosityVariable      = 1;        % Include effects of variable viscosity over space
bKappaVariable          = 1;        % Include effects of variable thermal conductivity

%% Output file name
name_file_out = ['output_data_' Test '_' num2str(num_grid_x) ...
    'x' num2str(num_grid_y) 'x' num2str(num_grid_z) '_KGP_Pt'];

%% Boundary conditions
BC_types = {'Periodic','Neumann','Dirichlet'}; % Select BC_Types{1},{2},{3}
% Define BC for BC_type & [u,v,w,P,T]
% Dirichlet: Set P<0 for impermeable boundary
% Dirichlet: Set T<=0 for no Temperature effect
BC.west   = {BC_types{1},0,0,0,0,0};
BC.east   = {BC_types{1},0,0,0,0,0};
BC.south  = {BC_types{1},0,0,0,0,0};
BC.north  = {BC_types{1},0,0,0,0,0};
BC.back   = {BC_types{1},0,0,0,0,0};
BC.front  = {BC_types{1},0,0,0,0,0};


%% 3D COMPRESSIBLE SOLVER
[u,v,w,E,rho,mu,kappa,c_v,c_p,P,T,ke,e,sos,Beta_T, Beta_v, Beta_s, Alpha_p,time,t_vec,ke_total,X,Y,Z, Invariants] = Compressible_Solver(Fluid,bSolver,x_0,y_0,z_0,L_x,L_y, L_z,...
        t_0,t_end,name_file_out,num_grid_x,num_grid_y,num_grid_z, CFL, n_iter_max, output_iter, RK_order, scheme, BC, ST, K, ...
        bViscous_5PointsStencil, bViscosityVariable, bKappaVariable, bEnergySplit, bEnergySplitEnthalpy, Test, Substance, HP_model, bPressureModel, bFlux, bFilter);

%% Data output
DataOutput(u,v,w,E,rho,mu,kappa,P,T,ke,e,c_p,c_v,sos,Beta_T,Beta_v,Beta_s,Alpha_p,X,Y,Z,t_vec,ke_total,Invariants,name_file_out);
    
  