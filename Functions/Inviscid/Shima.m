function F = Shima(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type, bSolver, Fluid, Substance)

% Preventing spurious pressure oscillations in split convective
% form discretization for compressible ﬂows
% JCP 2021
% Nao Shima, Yuichi Kuya, Yoshiharu Tamaki, Soshi Kawai∗

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
    K_L = 0.5*(u_L^2 + v_L^2 + w_L^2);
    e_L = E_L - K_L;
    K_R = 0.5*(u_R^2 + v_R^2 + w_R^2);
    e_R = E_R - K_R;
    F = 1/8*(rho_L + rho_R)*(u_L + u_R)*(K_L + K_R) + 1/4*(rho_L*e_L + rho_R*e_R)*(u_L + u_R) + 1/2*(u_L*P_R + u_R*P_L);

end



end