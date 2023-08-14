function [rhou_visc, rhov_visc, rhow_visc, rhoE_visc, Tau, q_div] = Calculate_Viscous(u,v,w,mu,k,T,dx,dy,dz,bViscosityVariable,bKappaVariable)  

    % Gradient computes de 2nd order derivative with respect to x,y,z
    [du_x,du_y,du_z]                = CentralDerivative_d1_2ndOrder(u);
    [dv_x,dv_y,dv_z]                = CentralDerivative_d1_2ndOrder(v);
    [dw_x,dw_y,dw_z]                = CentralDerivative_d1_2ndOrder(w);

    % Laplacian velocity
    [d2u_x,d2u_y,d2u_z]             = CentralDerivative_d2_2ndOrder(u);
    [d2v_x,d2v_y,d2v_z]             = CentralDerivative_d2_2ndOrder(v);
    [d2w_x,d2w_y,d2w_z]             = CentralDerivative_d2_2ndOrder(w);

    % Divergence velocity
    div_U = du_x./dx + dv_y./dy + dw_z./dz;
    
    % Derivative divergence
    [d_div_U_x,d_div_U_y,d_div_U_z]  = CentralDerivative_d1_2ndOrder(div_U);

    % Derivative mu
    [dmu_x,dmu_y,dmu_z]              = CentralDerivative_d1_2ndOrder(mu);

    if bViscosityVariable == 0
        dmu_x = 0.*dmu_x;
        dmu_y = 0.*dmu_y;
        dmu_z = 0.*dmu_z;
    end


    % Momentum equations
    rhou_visc = mu.*d2u_x./dx.^2 + mu.*d2u_y./dy.^2 + mu.*d2u_z./dz.^2 + ...
        + 1/3*mu.*d_div_U_x./dx + ...
        + dmu_x./dx.*(2*du_x./dx) + dmu_y./dy.*(du_y./dy + dv_x./dx) + dmu_z./dz.*(du_z./dz + dw_x./dx) + ...
        - 2/3*(dmu_x./dx).*(du_x./dx) - 2/3*(dmu_x./dx).*(dv_y./dy) - 2/3*(dmu_x./dx).*(dw_z./dz);

    rhov_visc = mu.*d2v_x./dx.^2 + mu.*d2v_y./dy.^2 + mu.*d2v_z./dz.^2 + ...
        + 1/3*mu.*d_div_U_y./dy + ...
        + dmu_x./dx.*(dv_x./dx + du_y./dy) + dmu_y./dy.*(2*dv_y./dy) + dmu_z./dz.*(dv_z./dz + dw_y./dy) + ...
        - 2/3*(dmu_y./dy).*(du_x./dx) - 2/3*(dmu_y./dy).*(dv_y./dy) - 2/3*(dmu_y./dy).*(dw_z./dz);

    rhow_visc = mu.*d2w_x./dx.^2 + mu.*d2w_y./dy.^2 + mu.*d2w_z./dz.^2 + ...
        + 1/3*mu.*d_div_U_z./dz + ...
        + dmu_x./dx.*(dw_x./dx + du_z./dz) + dmu_y./dy.*(dw_y./dy + dv_z./dz) + dmu_z./dz.*(2*dw_z./dz) + ...
        - 2/3*(dmu_z./dz).*(du_x./dx) - 2/3*(dmu_z./dz).*(dv_y./dz) - 2/3*(dmu_z./dz).*(dw_z./dz);

    % Calculate tau - Stress tensor
    Tau.Tau_xx = mu.*2.*du_x./dx - 2/3*mu.*div_U;
    Tau.Tau_xy = mu.*(du_y./dy + dv_x./dx);
    Tau.Tau_xz = mu.*(du_z./dz + dw_x./dx);
    Tau.Tau_yy = mu.*2.*dv_y./dy - 2/3*mu.*div_U;
    Tau.Tau_yx = mu.*(dv_x./dx + du_y./dy);
    Tau.Tau_yz = mu.*(dv_z./dz + dw_y./dy);
    Tau.Tau_zz = mu.*2.*dw_z./dz - 2/3*mu.*div_U;
    Tau.Tau_zx = mu.*(dw_x./dx + du_z./dz);
    Tau.Tau_zy = mu.*(dw_y./dy + dv_z./dz);

    % Fourier term (energy)
    [dT_x, dT_y, dT_z]                = CentralDerivative_d1_2ndOrder(T);
    [d2T_x,d2T_y,d2T_z]               = CentralDerivative_d2_2ndOrder(T);
    [dk_x, dk_y, dk_z]                = CentralDerivative_d1_2ndOrder(k);

    if bKappaVariable == 0
        dk_x = 0.*dk_x;
        dk_y = 0.*dk_y;
        dk_z = 0.*dk_z;
    end

    q_div = (dk_x./dx).*(dT_x./dx) + (dk_y./dy).*(dT_y./dy) + (dk_z./dz).*(dT_z./dz) + ...
        + k.*d2T_x./dx.^2 + k.*d2T_y./dy.^2 + k.*d2T_z./dz.^2;

    % Energy equation
    rhoE_visc = u.*rhou_visc + v.*rhov_visc + w.*rhow_visc + ...
        + (du_x./dx).*Tau.Tau_xx + (du_y./dy).*Tau.Tau_xy + (du_z./dz).*Tau.Tau_xz  + ...
        + (dv_x./dx).*Tau.Tau_yx + (dv_y./dy).*Tau.Tau_yy + (dv_z./dz).*Tau.Tau_yz  + ...
        + (dw_x./dx).*Tau.Tau_zx + (dw_y./dy).*Tau.Tau_zy + (dw_z./dz).*Tau.Tau_zz + ...
        + q_div;
    
    % Set outer points to 0
    rhou_visc([1,end],:,:)   = 0;
    rhou_visc(:,[1,end],:)   = 0;
    rhou_visc(:,:,[1,end])   = 0;
    rhov_visc([1,end],:,:)   = 0;
    rhov_visc(:,[1,end],:)   = 0;
    rhov_visc(:,:,[1,end])   = 0;
    rhow_visc([1,end],:,:)   = 0;
    rhow_visc(:,[1,end],:)   = 0;
    rhow_visc(:,:,[1,end])   = 0;
    rhoE_visc([1,end],:,:)   = 0;
    rhoE_visc(:,[1,end],:)   = 0;
    rhoE_visc(:,:,[1,end])   = 0;


end

