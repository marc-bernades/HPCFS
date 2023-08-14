function [u_exact,v_exact, rho_exact, P_exact] = Calculate_Exact_2D_Vortex(X,Y,Z,time,mu,rho,Fluid)

nu = mu./rho;

r         = sqrt((X - Fluid.Vortex_x0).^2 + (Y - Fluid.Vortex_y0).^2);
r_norm    = r/Fluid.r_v;
u_exact   = Fluid.U_0.*(1 - Fluid.Ma_vortex./Fluid.Ma*(((Y - Fluid.Vortex_y0)./Fluid.r_v).*exp((1 - r_norm.^2)/2)));
v_exact   = Fluid.U_0.*(Fluid.Ma_vortex./Fluid.Ma*(((X - Fluid.Vortex_x0)./Fluid.r_v).*exp((1 - r_norm.^2)/2)));
w_exact   = u_exact*0;
P_exact   = Fluid.P_0.*((1 - ((Fluid.gamma - 1)/2*Fluid.Ma_vortex^2).*exp((1 - r_norm.^2))).^(Fluid.gamma/(Fluid.gamma - 1)));
rho_exact = Fluid.rho_0*((1 - ((Fluid.gamma - 1)/2*Fluid.Ma_vortex^2).*exp((1 - r_norm.^2))).^(1/(Fluid.gamma - 1)));
T_exact   = P_exact./(rho_exact*Fluid.R_specific);


end

