function Flux = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance)

switch scheme
    case 'KGP'

        Flux = KGP(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type, bSolver, Fluid, Substance);

    case 'Shima'

        Flux = Shima(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type, bSolver, Fluid, Substance);

    case 'HLLC'

        Flux = HLLC(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type, bSolver, Fluid, Substance);
    
    case 'HLLC_plus'

        Flux = HLLC_plus(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type, bSolver, Fluid, Substance);

    otherwise
        disp('Please select inviscid scheme in flux form...')
        return


end



end