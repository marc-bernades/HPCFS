function F = KGP(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type, bSolver, Fluid, Substance)

% Kennedy, Gruber & Pirozzoli (KGP) scheme:
% G. Coppola , F. Capuano , S. Pirozzoli, L. de Luca.
% Numerically stable formulations of convective terms for turbulent compressible flows.
% Journal of Computational Physics, 382, 86-104, 2019.

F = 1/8*(rho_L + rho_R)*(u_L + u_R);

if var_type == 0
    F = F*2;

elseif var_type == 1
    F = F*(u_L + u_R) + 1/2*(P_L + P_R);

elseif var_type == 2
    F = F*(v_L + v_R);

elseif var_type == 3
    F = F*(w_L + w_R);

elseif var_type == 4
    F = F*(E_L + P_L/rho_L + E_R + P_R/rho_R);

end



end