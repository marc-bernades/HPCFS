function [S_L,S_R] = CalculateWavesSpeed(rho_L,rho_R,u_L,u_R,P_L,P_R,a_L,a_R, bSolver, Fluid, Substance)

% HLLC approximate Riemann solver:
% E. F. Toro.
% Riemann solvers and numerical methods for fluid dynamics.
% Springer, 2009.

% Obtain gamma
P_bar   = 0.5*(P_L + P_R);
rho_bar = 0.5*(rho_L + rho_R);
T_bar   = Calculate_T_fromPandRho( bSolver, P_bar,rho_bar,Fluid,Substance);
[~,~,gamma] = Calculate_SpecificHeatCapacities(bSolver, P_bar,rho_bar,T_bar, Fluid,Substance);
a_bar   = 0.5*( a_L + a_R );
P_pvrs  = 0.5*( P_L + P_R ) - 0.5*( u_R - u_L ).*rho_bar.*a_bar;
P_star  = max( 0.0, P_pvrs );
q_L     = 1.0;
if( P_star > P_L )
    q_L = sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_L ) - 1.0 ) );
end
q_R     = 1.0;
if( P_star > P_R )
    q_R = sqrt( 1.0 + ( ( gamma + 1.0 )/( 2.0*gamma ) )*( ( P_star/P_R ) - 1.0 ) );
end

% Wave speed estimates
S_L     = u_L - a_L*q_L;
S_R     = u_R + a_R*q_R;

% Direct wave speed estimates:
% B. Einfeldt.
% On Godunov-type methods for gas dynamics.
% SIAM Journal on Numerical Analysis, 25, 294-318, 1988.

hat_u = ( u_L*sqrt( rho_L ) + u_R*sqrt( rho_R ) )/( sqrt( rho_L ) + sqrt( rho_R ) );
hat_a = sqrt( ( ( a_L*a_L*sqrt( rho_L ) + a_R*a_R*sqrt( rho_R ) )/( sqrt( rho_L ) + sqrt( rho_R ) ) )+ 0.5*( ( sqrt( rho_L )*sqrt( rho_R ) )/( ( sqrt( rho_L ) + sqrt( rho_R ) )*( sqrt( rho_L ) + sqrt( rho_R ) ) ) )*( u_R - u_L )*( u_R - u_L ) );

S_L = min( u_L - a_L, hat_u - hat_a );
S_R = max( u_R + a_R, hat_u + hat_a );

end
