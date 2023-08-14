function [rhou_visc, rhov_visc, rhow_visc, rhoE_visc] = Calculate_Viscous_5PointStencil(u,v,w,mu,k,T,X,Y,Z)  

    % Delta x, y and z based on gradient
    [dx,~,~] = gradient(X);
    [~,dy,~] = gradient(Y);
    [~,~,dz] = gradient(Z);

    % Gradient computes de 2nd order derivative with respect to x,y,z
    [du_x,du_y,du_z]             = gradient(u);
    [dv_x,dv_y,dv_z]             = gradient(v);
    [dw_x,dw_y,dw_z]             = gradient(w);
    
    % Calculate tau - Stress tensor
    Tau_xx = mu.*2.*du_x./dx - 2/3*mu.*(du_x./dx + dv_y./dy + dw_z./dz);
    Tau_xy = mu.*(du_y./dy + dv_x./dx);
    Tau_xz = mu.*(du_z./dz + dw_x./dx);
    Tau_yy = mu.*2.*dv_y./dy - 2/3*mu.*(du_x./dx + dv_y./dy + dw_z./dz);
    Tau_yx = mu.*(dv_x./dx + du_y./dy);
    Tau_yz = mu.*(dv_z./dz + dw_y./dy);
    Tau_zz = mu.*2.*dw_z./dz - 2/3*mu.*(du_x./dx + dv_y./dy + dw_z./dz);
    Tau_zx = mu.*(dw_x./dx + du_z./dz);
    Tau_zy = mu.*(dw_y./dy + dv_z./dz);

    % Gradient tau
    [dTau_xx_x,dTau_xx_y,dTau_xx_z] = gradient(Tau_xx);
    [dTau_xy_x,dTau_xy_y,dTau_xy_z] = gradient(Tau_xy);
    [dTau_xz_x,dTau_xz_y,dTau_xz_z] = gradient(Tau_xz);
    [dTau_yy_x,dTau_yy_y,dTau_yy_z] = gradient(Tau_yy);
    [dTau_yx_x,dTau_yx_y,dTau_yx_z] = gradient(Tau_yx);
    [dTau_yz_x,dTau_yz_y,dTau_yz_z] = gradient(Tau_yz);
    [dTau_zz_x,dTau_zz_y,dTau_zz_z] = gradient(Tau_zz);
    [dTau_zx_x,dTau_zx_y,dTau_zx_z] = gradient(Tau_zx);
    [dTau_zy_x,dTau_zy_y,dTau_zy_z] = gradient(Tau_zy);

    % Viscous terms - Momentum
    rhou_visc = dTau_xx_x./dx + dTau_xy_y./dy + dTau_xz_z./dz;
    rhov_visc = dTau_yx_x./dx + dTau_yy_y./dy + dTau_yz_z./dz;
    rhow_visc = dTau_zx_x./dx + dTau_zy_y./dy + dTau_zz_z./dz;

    % Fourier term (energy)
    [dT_x,dT_y,dT_z]                = gradient(T);
    T_div                           = dT_x./dx + dT_y./dy + dT_z./dz;
    [dT_grad_x,dT_grad_y,dT_grad_z] = gradient(k.*T_div);
    q_div                           = dT_grad_x./dx + dT_grad_y./dy + dT_grad_z./dz;

    % Gradient tau*u_i
    [dTau_xx_u_x,dTau_xx_u_y,dTau_xx_u_z] = gradient(Tau_xx.*u);
    [dTau_xy_u_x,dTau_xy_u_y,dTau_xy_u_z] = gradient(Tau_xy.*u);
    [dTau_xz_u_x,dTau_xz_u_y,dTau_xz_u_z] = gradient(Tau_xz.*u);
    [dTau_yy_v_x,dTau_yy_v_y,dTau_yy_v_z] = gradient(Tau_yy.*v);
    [dTau_yx_v_x,dTau_yx_v_y,dTau_yx_v_z] = gradient(Tau_yx.*v);
    [dTau_yz_v_x,dTau_yz_v_y,dTau_yz_v_z] = gradient(Tau_yz.*v);
    [dTau_zz_w_x,dTau_zz_w_y,dTau_zz_w_z] = gradient(Tau_zz.*w);
    [dTau_zx_w_x,dTau_zx_w_y,dTau_zx_w_z] = gradient(Tau_zx.*w);
    [dTau_zy_w_x,dTau_zy_w_y,dTau_zy_w_z] = gradient(Tau_zy.*w);

    % Viscous term - Energy
    rhoE_visc = dTau_xx_u_x./dx + dTau_xy_u_y./dy + dTau_xz_u_z./dz + ...
        dTau_yx_v_x./dx + dTau_yy_v_y./dy + dTau_yz_v_z./dz + ...
         dTau_zx_w_x./dx + dTau_zy_w_y./dy + dTau_zz_w_z./dz + ...
         q_div;
    rhoE_visc = q_div;

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

