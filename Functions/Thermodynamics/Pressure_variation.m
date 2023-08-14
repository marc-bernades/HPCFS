function dP = Pressure_variation(P,rho,u,v,w,sos,dx,dy,dz,c_v,Alpha_p,Beta_T,Tau,q_div)

% Find P from Kawai 2015 implementation non-conservative general pressure
% equation (Eq.19)
% Divergence u
[du_x,du_y,du_z] = CentralDerivative_d1_2ndOrder(u);
[dv_x,dv_y,dv_z] = CentralDerivative_d1_2ndOrder(v);
[dw_x,dw_y,dw_z] = CentralDerivative_d1_2ndOrder(w);

divU = du_x./(dx) + dv_y./(dy) + dw_z./(dz);


% Gradient P
[dP_x,dP_y,dP_z] = CentralDerivative_d1_2ndOrder(P);

u_GradP = u.*(dP_x./dx) + v.*(dP_y./dy) + w.*(dP_z./dz);

% First we solve current dP based on Treshima 2012 Equation 18 > Ammended
% Kawai 2015 Eq 19
% Inviscid case ONLY
dP_inv = -u_GradP - rho.*sos.^2.*divU;

% Conservative case for IDEAL GAS GAMMA = 1
% Pu_conv_CD   = Inviscid_CD(P,u,v,w,ones(size(rho)),dx,dy,dz);
% dP_inv       = -Pu_conv_CD;

% Viscous part
Tau_dU = (du_x./dx).*Tau.Tau_xx + (du_y./dy).*Tau.Tau_xy + (du_z./dz).*Tau.Tau_xz  + ...
        + (dv_x./dx).*Tau.Tau_yx + (dv_y./dy).*Tau.Tau_yy + (dv_z./dz).*Tau.Tau_yz  + ...
        + (dw_x./dx).*Tau.Tau_zx + (dw_y./dy).*Tau.Tau_zy + (dw_z./dz).*Tau.Tau_zz;
dP_visc = Alpha_p./(c_v.*Beta_T).*(1./rho.*(Tau_dU - q_div));

% Final
dP = dP_inv + dP_visc;

%% Equation based on ideal gas
% Gradient P
% [dPu_x,dPu_y,dPu_z] = CentralDerivative_d1_2ndOrder(P.*u);
% [dPv_x,dPv_y,dPv_z] = CentralDerivative_d1_2ndOrder(P.*v);
% [dPw_x,dPw_y,dPw_z] = CentralDerivative_d1_2ndOrder(P.*w);
% divPU = dPu_x./(dx) + dPv_y./(dy) + dPw_z./(dz);


% Pu_conv_CD = Inviscid_CD(P,u,v,w,ones(size(rho)),dx,dy,dz);
% Pu_conv_Cu = Inviscid_Cu(P,u,v,w,ones(size(rho)),dx,dy,dz);
% divPU      = 0.5*Pu_conv_CD + 0.5*Pu_conv_Cu;

% Gamma = 1.4;
% dP = - divPU - (Gamma - 1).*P.*divU;
% 
% 
% 
% [dP_x,dP_y,dP_z] = CentralDerivative_d1_2ndOrder(P);
% 
% uGradP = u.*dP_x./(dx) + v.*dP_y./(dy) + w.*dP_z./(dz);
% dP = -uGradP - Gamma*P.*divU;

% Set outer points to 0
dP([1,end],:,:)    = 0;
dP(:,[1,end],:)    = 0;
dP(:,:,[1,end])    = 0;


end

