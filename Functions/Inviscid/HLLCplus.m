function Flux = HLLCplus(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type)
% HLLC-type Riemann solver for all-speed flows:
% S. Chen, B. Lin, Y. Li, C. Yan.
% HLLC+: low-Mach shock-stable HLLC-type Riemann solver for all-speed flows.
% SIAM Journal of Scientific Computing, 4, B921-B950, 2020.

F_L = rho_L.*u_L; F_R = rho_R.*u_R;
U_L = rho_L;      U_R = rho_R;

if var_type == 0
    F_L = F_L; FR = F_R;
    U_L = U_L; UR = U_R;

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
[S_L,S_R] = CalculateWavesSpeed(rho_L,rho_R,u_L,u_R,P_L,P_R,a_L,a_R);

phi_L = rho_L./(S_L - u_L);
phi_R = rho_R./(S_R - u_R);
S_star = (P_R - P_L + phi_L.*u_L - phi_R.*u_R)./(phi_L - phi_R);
U_star_L = rho_L.*(S_L - u_L)./(S_L - S_star);
U_star_R = rho_R.*(S_R - u_R)./(S_R - S_star);
M = min(A1,max(1./a_L.*sqrt(u_L.^2 + v_L.^2 + w_L.^2),1./a_R.*sqrt(u_R.^2 + v_R.^2 + w_R.^2)));
f_M = M.*sqrt(4 + (1- M.^2).^2)./(1 + M.^2);
h = min(P_L./P_R, P_R./P_L);
g = 1 - h.^M;
A_p_L = phi_L.*phi_R./(phi_R - phi_L);
A_p_R = phi_L.*phi_R./(phi_R - phi_L);

if var_type == 0
    U_star_L = U_star_L; U_star_R = U_star_R;
    A_p_L    = 0;        A_p_R    = 0;

elseif var_type == 1
    U_star_L = U_star_L.*S_star;              U_star_R = U_star_R.*S_star;
    A_p_L    = A_p_L.*(f_M - 1).*(u_R - u_L); A_p_R    = A_p_R.*(f_M - 1).*(u_R - u_L);   

elseif var_type == 2
    U_star_L = U_star_L.*v_L;              
    U_star_R = U_star_R.*v_R;
    A_p_L    = A_p_L.*S_L./(S_L - S_star).*g.*(v_R - v_L);
    A_p_R    = A_p_R.*S_R./(S_R - S_star).*g.*(v_R - v_L);

elseif var_type == 3
    U_star_L = U_star_L.*w_L;              
    U_star_R = U_star_R.*w_R;
    A_p_L    = A_p_L.*S_L./(S_L - S_star).*g.*(w_R - w_L);
    A_p_R    = A_p_R.*S_R./(S_R - S_star).*g.*(w_R - w_L);

elseif var_type == 4
    U_star_L = U_star_L.*(E_L + (S_star - u_L).*(S_star + P_L./(rho_L.*(S_L - u_L))));              
    U_star_R = U_star_R.*(E_R + (S_star - u_R).*(S_star + P_R./(rho_R.*(S_R - u_R))));
    A_p_L    = A_p_L.*(f_M - 1).*(u_R - u_L).*S_star;
    A_p_R    = A_p_R.*(f_M - 1).*(u_R - u_L).*S_star;

end

F_star_L = F_L + S_L.*(U_star_L - U_L) + A_p_L;
F_star_R = F_R + S_R.*(U_star_R - U_R) + A_p_R;

F = 0;
if 0 <= S_L
    F = F_L;
elseif S_L <= 0 && 0 <= S_star
    F = F_star_L;
elseif S_star <= 0 && 0 <= S_R
    F = F_star_R;
elseif 0 >= S_R
    F = F_R;
end





end