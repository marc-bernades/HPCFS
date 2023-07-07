function  dt = Calculate_TimeStep(bSolver,rho,rhou,rhov,rhow,CFL,sos,gamma,c_p,kappa,mu,num_grid,L,dx,dy,dz,Fluid,Test)

if strcmp(bSolver,'Real')

    u = rhou./rho;
    v = rhov./rho;
    w = rhow./rho;

    % Guess Pr number
    Pr       = gamma;

    % Calculation Pr number
    epsilon  = 10E-15;
    Cond     = kappa > epsilon;
    Pr(Cond) = c_p(Cond).*mu(Cond)./kappa(Cond);

    % x-direction inviscid and viscous terms
    if sum(u) ~= 0
        S_x = abs(u) + sos;
        dt = min(min(min(CFL*dx./S_x)));
        dt = min(dt,min(min(min(CFL.*Pr.*rho.*dx.^2./(mu.*gamma + epsilon)))));
    end

    if sum(v) ~= 0
        % y-direction inviscid and viscous terms
        S_y = abs(v) + sos;
        dt = min(dt,min(min(min(CFL*dy./S_y))));
        dt = min(dt,min(min(min(CFL.*Pr.*rho.*dy.^2./(mu.*gamma + epsilon)))));
    end

    if sum(w) ~= 0
        % z-direction inviscid and viscous terms
        S_z = abs(w) + sos;
        dt = min(dt,min(min(min(CFL*dz./S_z))));
        dt = min(dt,min(min(min(CFL.*Pr.*rho.*dz.^2./(mu.*gamma + epsilon)))));
    end

    if sum(u) == 0 & sum(v) == 0 & sum(w) == 0
        u_vec = sqrt((rhou.*rhou + rhov.*rhov + rhow.*rhow)./rho.^2);
        U_ref = max(max(max(max(max(u_vec))),max(max(max(u_vec + sos)))),max(max(max(u_vec - sos))));
        dt    = CFL*(L/num_grid)/U_ref;
    end
else

    if strcmp(Test,'1D_Adv')
        U_ref = Fluid.U_0;
    else
        U_ref = sqrt(max(max(max((rhou.*rhou + rhov.*rhov + rhow.*rhow)./rho.^2)))) + sqrt(Fluid.gamma*Fluid.R_specific*Fluid.T_0);
    end
    %U_ref = U_0 + sqrt(gamma*R_specific*T_0);
    dt    = CFL*(L/num_grid)/U_ref;
    %dt    = 1e-4;
    %dt/((2*pi/num_grid)/U_ref)

end
end

