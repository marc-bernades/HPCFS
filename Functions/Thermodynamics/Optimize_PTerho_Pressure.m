function [P,T] = Optimize_PTerho_Pressure(P_current,rho_target,Substance,u,v,w,dx,dy,dz,dt,...
    c_v,sos,Beta_T, Alpha_p,Tau,q_div)

% Divergence u
[du_x,du_y,du_z] = CentralDerivative_d1_2ndOrder_mex(u);
[dv_x,dv_y,dv_z] = CentralDerivative_d1_2ndOrder_mex(v);
[dw_x,dw_y,dw_z] = CentralDerivative_d1_2ndOrder_mex(w);

divU = du_x./(dx) + dv_y./(dy) + dw_z./(dz);


% Gradient P
[dP_x,dP_y,dP_z] = CentralDerivative_d1_2ndOrder_mex(P_current);

u_GradP = u.*(dP_x./dx + dP_y./dy + dP_z./dz) + ...
    v.*(dP_x./dx + dP_y./dy + dP_z./dz) + ...
    w.*(dP_x./dx + dP_y./dy + dP_z./dz);

% Viscous part of internal energy (Tau : nabla(u'))
Tau_GradU = Tau.Tau_xx.*du_x./(dx) + Tau.Tau_xy.*du_y./(dy) + Tau.Tau_xz.*du_z./(dz) + ...
    Tau.Tau_yx.*dv_x./(dx) + Tau.Tau_yy.*dv_y./(dy) + Tau.Tau_yz.*dv_z./(dz) + ...
    Tau.Tau_zx.*dw_x./(dx) + Tau.Tau_zy.*dw_y./(dy) + Tau.Tau_zz.*dw_z./(dz);

% First we solve current dP based on Treshima 2012 Equation 18 > Ammended
% Kawai 2015 Eq 19
dP = -u_GradP - rho_target.*sos.^2.*divU + Alpha_p./(c_v.*Beta_T).*(1./rho_target).*(Tau_GradU + q_div);

% Euler method considering explicit depending on current time
P  = P_current + dt*dP;

% We can update temperature based on P and rho
T = Calculate_T_fromPandRho( P,rho_target,Substance);

end
