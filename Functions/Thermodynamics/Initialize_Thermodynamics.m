function [rho,e,ke,E,sos,c_v,c_p,gamma,Beta_T, Beta_v, Beta_s, Alpha_p] = Initialize_Thermodynamics(bSolver,P,T,u,v,w,bPressureModel,Fluid,Substance)

rho             = Calculate_Rho_from_TandP(bSolver, T, P, Fluid, Substance);
e               = Calculate_e_from_TandPandRho(bSolver, T,rho,P, Fluid, Substance );
ke              = 0.5*(u.^2 + v.^2 + w.^2);
E               = e + ke;

[c_p,c_v,gamma]      = Calculate_SpecificHeatCapacities(bSolver, P,rho,T, Fluid, Substance);
[sos, Beta_T, Beta_v, Beta_s, Alpha_p]  = Calculate_sos(bSolver,rho,T,c_p,P,bPressureModel,Fluid,Substance);

end




