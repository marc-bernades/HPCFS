%% HIGH PRESSURE - RATIO DENSITY SCALE OVER KOLMOGOROV SCALE
% Load results
DATA        = readtable('output_data_1D_HighPressure_50x1x1_2Pc_0.csv');

% Index only unique Y and Z
Index_Z = unique(DATA.Z); Index_Z = find(DATA.Z == Index_Z(2));
Index   = Index_Z(5:3:end-3);

% Variables
T     = DATA.T(Index);
rho   = DATA.rho(Index);
c_p   = DATA.c_p(Index);
mu    = DATA.mu(Index);
kappa = DATA.kappa(Index);

% Calculate volume expansivity
% PengRobinson
Substance = 'N2';
[a,b,R,dadT,d2adT2,NASA_coefficients] = PengRobinson(T, Substance);
% Volume
v = 1./rho;
% Thermodynamic expressions for sos
dP_dv_constant_T = -1*R*T./(v - b).^2 + a.*(2*v + 2*b)./(v.^2 + 2*v*b - b^2).^2;
dP_dT_constant_v = R./(v - b) - dadT./(v.^2 + 2.*v.*b - b.^2);
Beta_T           = -1./(v.*dP_dv_constant_T); %Isothermal compressability
Beta_v           = -1*(dP_dT_constant_v./(v.*dP_dv_constant_T)); %Volume expansivity

% Temperature variation across Pb
Beta_v_pb = max(Beta_v);
Beta_v_int = cumtrapz(T,Beta_v); %Integral from T1 to T2
deltaT = 1/Beta_v_pb*Beta_v_int(end);

% Prandtl number at Pb
Pr    = c_p.*mu./kappa; 
Pr_pb = max(Pr);

% Reynolds number of the large-scale relative motion in mixing layer
% Define channel flow height
delta = 100*1E-6;
H     = 2*delta;
% Define bulk values
rho_b = 502.03;
u_b   = 2.434;
mu_b  = 3.972299*1E-5;

Re_M = rho_b*u_b*H/mu_b;

% Ratio density gradient over Kolmogorov scales
% Near pb
Ratio = deltaT*Re_M^(1/4)/((T(end) - T(1))*Pr_pb^0.5);
% Mixing layer
Ratio_2 = ((deltaT/(T(end) - T(1)))^3)*Re_M^(3/4);