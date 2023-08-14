function [rhou, rhov, rhow, rhoE] = Update_ConservedVariables(rho,u,v,w,E)
    rhou = rho.*u;
    rhov = rho.*v;
    rhow = rho.*w;
    rhoE = rho.*E;
end

