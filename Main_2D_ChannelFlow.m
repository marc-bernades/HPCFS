%% 3D Compressible solver %%
close all; clear all; clc

%% SET PROBLEM PARAMETERS %%
% Thermodynamics framework
bSolver           = 'Real';                    % Define Ideal or Real Gas Model
% Substance
Substance         = 'N2';                      % Fluid selected substance
HP_model          = 'HighPressure';                % Transport coefficients model: 'Constant', 'LowPressure', 'HighPressure'

% Fluid properties
if strcmp(bSolver,'Ideal')
    % IDEAL

else
    % REAL
    Fluid.delta     = 100*1E-6; % Half channel height
    [~, T_c, P_c, ~, ~, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Substance);
    Fluid.P_0       = 2*P_c;                    % Rereference Pressure [Pa]
    Fluid.T_bw      = 0.75*T_c;
    Fluid.T_tw      = 1.5*T_c;
    Fluid.Re_tau_bw = 100;
    Fluid.rho_bw    = Calculate_Rho_from_TandP( bSolver, Fluid.T_bw,Fluid.P_0, Fluid, Substance ); %839.39;
    Fluid.mu_bw     = Calculate_HighPressure_Transport_Coeff(bSolver, 0, 0, Fluid, Substance, Fluid.T_bw, Fluid.rho_bw, HP_model); % 0.00016;
    Fluid.u_tau_bw  = Fluid.mu_bw*Fluid.Re_tau_bw/(Fluid.rho_bw*Fluid.delta);
    Fluid.tau_bw    = Fluid.rho_bw*Fluid.u_tau_bw^2;
    
    % Perturbations initialization
    Fluid.kappa_vK  = 0.41; % Von Karman constant
    Fluid.y_0       = Fluid.mu_bw/Fluid.rho_bw/(9*Fluid.u_tau_bw); % Smooth-wall roughness bottom wall [m]
    Fluid.u_0       = ( Fluid.u_tau_bw/Fluid.kappa_vK )*( log( Fluid.delta/Fluid.y_0 ) + ( Fluid.y_0/Fluid.delta ) - 1.0 ); % Volume-average of a log-law velocity profile [m/s]
    Fluid.alpha     = 0.25; % Magnitude perturbations


end

% Problem parameters
x_0           = 0;                          % Domain origin in x-direction [m]
y_0           = 0;                          % Domain origin in y-direction [m]
z_0           = 0;                          % Domain origin in z-direction [m]
L_x           = 4*pi*Fluid.delta;           % Size of domain in x-direction [m]
L_y           = Fluid.delta*2;              % Size of domain in y-direction [m]
L_z           = 0.01*1;                     % Size of domain in z-direction [m]
t_0           = 0.0;                        % Initial time [s]
t_c           = Fluid.delta/Fluid.u_tau_bw; % Reference time [s]
t_end         = 1*t_c;                      % Final time [s]
Test          = '2D_ChannelFlow';         

% Computational parameters
num_grid_x  = 32;
num_grid_y  = 64;
num_grid_z  = 1;
CFL         = 0.9;                  % CFL time-step stability coefficient
n_iter_max  = 1E10;                 % Maximum number of time iterations
output_iter = 0:t_c/100:t_end;      % Output data every selected time step
RK_order    = 4;                    % Time-integration Runge-Kutta order

% Stretching
% Factors: x = x_0 + L*eta + A*( 0.5*L - L*eta )*( 1.0 - eta )*eta, with eta = ( l - 0.5 )/num_grid  A < 0: stretching at ends; A = 0: uniform; A > 0: stretching at center
% Calculate factors
y_plus = 0.1;
Re_tau = Fluid.Re_tau_bw;
delta  = 1; % Unitary for A calculation
u_tau  = 1;

Fluid.A_x = 0;
Fluid.A_y = Calculate_stretching_factor(y_plus, u_tau, Re_tau, delta,num_grid_y, y_0, 2*delta); %-1.9166884139482563;
Fluid.A_z = 0;

% Pressure spurious oscillation
bPressureModel          = 1;        % Use pressure model to avoid spurious oscillations

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
scheme.e        = [1/4 1/4 1/4 1/4 0]; % Shima [1/2 0 1/2 0 0];
scheme.k        = [1/4 1/4 1/4 1/4 0];
scheme.p        = [1/4 1/4 1/4 1/4 0]; % KGP [1/2 0 1/2 0 0], KGP with enthalpy [1/4 1/4 1/4 1/4 0], Shima [0 0 1 0 0];  
  
% Viscous Scheme
bViscous_5PointsStencil = 0;        % Gate for 5-points stencil choice
bViscosityVariable      = 1;        % Include effects of variable viscosity over space
bKappaVariable          = 1;        % Include effects of variable thermal conductivity

% Output file name
name_file_out = ['output_data_' Test '_' num2str(num_grid_x) ...
    'x' num2str(num_grid_y) 'x' num2str(num_grid_z) '_KGP_Pt'];

% Boundary conditions
BC_types = {'Periodic','Neumann','Dirichlet'}; % Select BC_Types{1},{2},{3}
% Define BC for BC_type & [u,v,w,P,T]
% Dirichlet: Set u,v,w=-1 for Symmetry and <-1 for Neuman and -1<u<0 for
% periodic
% Dirichlet: Set P<0 for impermeable boundary
% Dirichlet: Set T<=0 for impermeable boundary
BC.west   = {BC_types{1},0,0,0,0,0};
BC.east   = {BC_types{1},0,0,0,0,0};
BC.south  = {BC_types{3},0,0,0,-1,0.75*T_c};
BC.north  = {BC_types{3},0,0,0,-1,1.5*T_c};
BC.back   = {BC_types{1},0,0,0,0,0};
BC.front  = {BC_types{1},0,0,0,0,0};

% Source terms
ST.f_rhou = 0;
ST.f_rhov = 0;
ST.f_rhow = 0;
ST.f_rhoE = 0;

% Define Controller
K.b_active      = 1;    % Boolean activating controller
K.case          = Test; % Define case to apply controller algorithm
K.kp            = 0.1;  % Proportional
K.tau_bw_target = Fluid.tau_bw; % Target
K.mu_bw         = Fluid.mu_bw;
K.controller_output = K.tau_bw_target/Fluid.delta; %Initial guess
ST.f_rhou           = K.controller_output; % Set initial Source term value


%% 3D COMPRESSIBLE SOLVER
[u,v,w,E,rho,mu,kappa,c_v,c_p,P,T,ke,e,sos,Beta_T, Beta_v, Beta_s, Alpha_p,time,t_vec,ke_total,X,Y,Z, Invariants] = Compressible_Solver(Fluid,bSolver,x_0,y_0,z_0,L_x,L_y, L_z,...
        t_0,t_end,name_file_out,num_grid_x,num_grid_y,num_grid_z, CFL, n_iter_max, output_iter, RK_order, scheme, BC, ST, K, ...
        bViscous_5PointsStencil, bViscosityVariable, bKappaVariable, bEnergySplit, bEnergySplitEnthalpy, Test, Substance, HP_model, bPressureModel, bFlux, bFilter);

%% Data output
DataOutput(u,v,w,E,rho,mu,kappa,P,T,ke,e,c_p,c_v,sos,Beta_T,Beta_v,Beta_s,Alpha_p,X,Y,Z,t_vec,ke_total,Invariants,name_file_out);
    
  