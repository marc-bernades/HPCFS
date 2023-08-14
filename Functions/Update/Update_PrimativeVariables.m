function [u, v, w, E] = Update_PrimativeVariables(rho,rhou,rhov,rhow,rhoE)
    
    u = rhou./rho;
    v = rhov./rho;
    w = rhow./rho;
    E = rhoE./rho;

end

