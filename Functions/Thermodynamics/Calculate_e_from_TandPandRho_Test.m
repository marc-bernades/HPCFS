function [ e ] = Calculate_e_from_TandPandRho_Test(bSolver, T,rho,P,Fluid, Substance,c_p,Beta_v )

if strcmp(bSolver,'Real')
    e = T.*(c_p.*rho - Beta_v.*P)./rho; 


else % IDEAL GAS
    e     = P./(rho*(Fluid.gamma - 1));


end
end

