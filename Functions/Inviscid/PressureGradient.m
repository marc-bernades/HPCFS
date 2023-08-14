function [dP_rhou, dP_rhov, dP_rhow, dP_rhoE] = PressureGradient(u,v,w,P,dx,dy,dz)

    % CentralDerivative_d1_2ndOrder computes de 2nd order derivative with respect to x,y,z
    [dP_x,dP_y,dP_z]        = CentralDerivative_d1_2ndOrder(P);
    dP_rhou = dP_x./(dx);
    dP_rhov = dP_y./(dy);
    dP_rhow = dP_z./(dz);

    [dPu_x,~,~]             = CentralDerivative_d1_2ndOrder(P.*u);
    [~,dPv_y,~]             = CentralDerivative_d1_2ndOrder(P.*v);
    [~,~,dPw_z]             = CentralDerivative_d1_2ndOrder(P.*w);

    dP_rhoE = dPu_x./(dx) + dPv_y./(dy) + dPw_z./(dz);

end

