function F = HLLC(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type, bSolver, Fluid, Substance)

% Harten-Lax-van Leer-Contact (HLLC) Riemman solver:
% E. F. Toro, M. Spruce, W. Speares.
% Restoration of the contact surface in the HLL-Riemann solver.
% Shock Waves, 4, 25-34, 1994.

F_L = rho_L.*u_L; F_R = rho_R.*u_R;
U_L = rho_L;      U_R = rho_R;

if var_type == 0
    F_L = F_L*1.0; F_R = F_R*1.0;
    U_L = U_L*1.0; U_R = U_R*1.0;

elseif var_type == 1
    F_L = F_L.*u_L + P_L; F_R = F_R.*u_R + P_R;
    U_L = U_L.*u_L;       U_R = U_R.*u_R;

elseif var_type == 2
    F_L = F_L.*v_L;       F_R = F_R.*v_R;
    U_L = U_L.*v_L;       U_R = U_R.*v_R;

elseif var_type == 3
    F_L = F_L.*w_L;       F_R = F_R.*w_R;
    U_L = U_L.*w_L;       U_R = U_R.*w_R;

elseif var_type == 4
    F_L = F_L.*E_L + u_L.*P_L;  F_R = F_R.*E_R + u_R.*P_R;
    U_L = U_L.*E_L;             U_R = U_R.*E_R;

end

% CalculateWavesSpeed
[S_L,S_R] = CalculateWavesSpeed(rho_L,rho_R,u_L,u_R,P_L,P_R,a_L,a_R, bSolver, Fluid, Substance);


S_star = (P_R - P_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R))./(rho_L*(S_L - u_L) - rho_R*(S_R - u_R));
U_star_L = rho_L.*(S_L - u_L)./(S_L - S_star);
U_star_R = rho_R.*(S_R - u_R)./(S_R - S_star);


if var_type == 0
    U_star_L = U_star_L*1.0; U_star_R = U_star_R*1.0;

elseif var_type == 1
    U_star_L = U_star_L.*S_star;
    U_star_R = U_star_R.*S_star;

elseif var_type == 2
    U_star_L = U_star_L.*v_L;
    U_star_R = U_star_R.*v_R;

elseif var_type == 3
    U_star_L = U_star_L.*w_L;
    U_star_R = U_star_R.*w_R;

elseif var_type == 4
    U_star_L = U_star_L.*(E_L + (S_star - u_L).*(S_star + P_L./(rho_L.*(S_L - u_L))));
    U_star_R = U_star_R.*(E_R + (S_star - u_R).*(S_star + P_R./(rho_R.*(S_R - u_R))));

end

F_star_L = F_L + S_L.*(U_star_L - U_L);
F_star_R = F_R + S_R.*(U_star_R - U_R);

F = 0;
if 0 <= S_L
    F = F_L;
elseif (S_L <= 0) && (0 <= S_star)
    F = F_star_L;
elseif (S_star <= 0) && (0 <= S_R)
    F = F_star_R;
elseif 0 >= S_R
    F = F_R;
end




end