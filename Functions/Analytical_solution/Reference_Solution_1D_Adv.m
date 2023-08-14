function [rho_ref, P_ref, T_ref, e_ref, u_ref, v_ref, w_ref, ke_ref, E_ref ] = Reference_Solution_1D_Adv(X,Fluid,Substance)
% Load Fluid constants
rho_min = Fluid.rho_min;
rho_max = Fluid.rho_max;
P_0     = Fluid.P_0;
U_0     = Fluid.U_0;
Case    = Fluid.Case;

% Initial rho based on Ma et al. 2017 JCP
% rho_ref            = zeros(size(X));
P_ref            = zeros(size(X));

switch Case
    case 'Sharp'
        Cond_max       = (0.25 < X) & (X < 0.75);
        rho_ref(Cond_max)  = rho_max;
        rho_ref(~Cond_max) = rho_min;
    case 'Smooth'
        rho_ref            = (rho_min + rho_max)/2 + (rho_max - rho_min)/2*sin(2*pi*(X));
%         Cond_max = (X > 0.25) & (X <= 0.75);
%         rho_ref(Cond_max)            = (rho_min + rho_max)/2 + (rho_max - rho_min)/2*sin(2*pi*(X(Cond_max) - 0.25));
%         rho_ref(~Cond_max)           = (rho_min + rho_max)/2;
%         P_ref(Cond_max)            = P_0 + P_0/2*sin(2*pi*(X(Cond_max) - 0.25));
%         P_ref(~Cond_max)           = P_0;
%         rho_ref            = (rho_min + rho_max)/2 - (rho_max - rho_min)/2*tanh(6*pi*(X - 0.5));
%         rho_ref            = (rho_min + rho_max)/2 + (rho_max - rho_min)/2*sin(2*pi*X);
%         rho_ref            = (rho_min + rho_max)/2 + (rho_max - rho_min)/2*sin(6*pi*(X - 0.5));

end

% Constant P based on P_0
P_ref   = P_0 + zeros(size(X));
% rho_ref   = rho_max + zeros(size(X)); 

% Constant u based on U_0
u_ref   = U_0 + zeros(size(X));
v_ref   = zeros(size(X));
w_ref   = zeros(size(X));

% Iterate T within EOS
T_ref   = Calculate_T_fromPandRho( 'Ideal',P_ref,rho_ref,Fluid,Substance);

% Calculate e based on T, P and rho
e_ref   = Calculate_e_from_TandPandRho(  'Ideal',T_ref,rho_ref,P_ref,Fluid,Substance );
ke_ref  = 0.5*u_ref.^2;
E_ref   = e_ref + ke_ref;


end