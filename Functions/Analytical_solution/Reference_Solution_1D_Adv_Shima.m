function [rho_ref, P_ref, T_ref, e_ref, u_ref, v_ref, w_ref, ke_ref, E_ref ] = Reference_Solution_1D_Adv_Shima(bSolver,X,Fluid,Substance)
% Constant P based on P_0
P_0     = Fluid.P_0;
P_ref   = P_0 + zeros(size(X)); 

% Constant u based on U_0
U_0     = Fluid.U_0;
u_ref   = U_0 + zeros(size(X));
v_ref   = zeros(size(X));
w_ref   = zeros(size(X));

% Final time
if isfield(Fluid,'t_final')
    t = Fluid.t_final;
else
    t = 0;
end

% Smooth part > Non-linear rho or T input
if isfield(Fluid,'rho_min')
    rho_min = Fluid.rho_min;
    rho_max = Fluid.rho_max;
    % Smooth test proposed by Shima 2021
    rho_ref = (rho_min + rho_max)/2  + (rho_max - rho_min)/2*sin(2*pi*(X - U_0*t)); % Correct for final advected time
    % Iterate T within EOS
    T_ref   = Calculate_T_fromPandRho( bSolver,P_ref,rho_ref,Fluid,Substance);
else
    T_min = Fluid.T0_min;
    T_max = Fluid.T0_max;
    T_ref = (T_min + T_max)/2  + (T_max - T_min)/2*sin(2*pi*X);
    rho_ref   = Calculate_Rho_from_TandP( bSolver,T_ref,P_ref,Fluid,Substance);

end





% Calculate e based on T, P and rho
e_ref   = Calculate_e_from_TandPandRho( bSolver,T_ref,rho_ref,P_ref, Fluid, Substance );
ke_ref  = 0.5*u_ref.^2;
E_ref   = e_ref + ke_ref;


end