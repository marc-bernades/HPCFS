function [C_D, C_u, C_fi, C_rho, C_L] = ConvectiveTerm(rho,u,v,w,E,rhou,rhov,rhow,rhoE,dx,dy,dz)
    

    % CentralDerivative_d1_2ndOrder computes de 2nd order derivative with respect to x,y,z
    [du_x,du_y,du_z]             = CentralDerivative_d1_2ndOrder(u);
    [dv_x,dv_y,dv_z]             = CentralDerivative_d1_2ndOrder(v);
    [dw_x,dw_y,dw_z]             = CentralDerivative_d1_2ndOrder(w);
    [drho_x,drho_y,drho_z]       = CentralDerivative_d1_2ndOrder(rho);
    [dE_x,  dE_y,  dE_z]         = CentralDerivative_d1_2ndOrder(E);

    [duu_x,~,~]                  = CentralDerivative_d1_2ndOrder(u.*u);
    [~,    dvv_y,~]              = CentralDerivative_d1_2ndOrder(v.*v);
    [~,    ~,    dww_z]          = CentralDerivative_d1_2ndOrder(w.*w);

    [duv_x,duv_y,~]              = CentralDerivative_d1_2ndOrder(u.*v);
    [duw_x,~,    duw_z]          = CentralDerivative_d1_2ndOrder(u.*w);
    [~,    dvw_y,dvw_z]          = CentralDerivative_d1_2ndOrder(v.*w);

    [duE_x,~,    ~]              = CentralDerivative_d1_2ndOrder(u.*E);
    [~,    dvE_y,~]              = CentralDerivative_d1_2ndOrder(v.*E);
    [~,    ~,    dwE_z]          = CentralDerivative_d1_2ndOrder(w.*E);

    [drhou_x,drhou_y,drhou_z]    = CentralDerivative_d1_2ndOrder(rhou);
    [drhov_x,drhov_y,drhov_z]    = CentralDerivative_d1_2ndOrder(rhov);
    [drhow_x,drhow_y,drhow_z]    = CentralDerivative_d1_2ndOrder(rhow);
    [drhoE_x,drhoE_y,drhoE_z]    = CentralDerivative_d1_2ndOrder(rhoE);

    [drhouu_x,~,       ~]        = CentralDerivative_d1_2ndOrder(rhou.*u);
    [~,       drhovv_y,~]        = CentralDerivative_d1_2ndOrder(rhov.*v);
    [~,        ~,      drhoww_z] = CentralDerivative_d1_2ndOrder(rhow.*w);

    [drhouv_x,drhouv_y,~]        = CentralDerivative_d1_2ndOrder(rhou.*v);
    [drhouw_x,~,       drhouw_z] = CentralDerivative_d1_2ndOrder(rhou.*w);
    [~,       drhovw_y,drhovw_z] = CentralDerivative_d1_2ndOrder(rhov.*w);

    [drhouE_x,~       ,~]        = CentralDerivative_d1_2ndOrder(rhou.*E);
    [~       ,drhovE_y,~]        = CentralDerivative_d1_2ndOrder(rhov.*E);
    [~       ,~       ,drhowE_z] = CentralDerivative_d1_2ndOrder(rhow.*E);



    %% Convective terms
    % CD - Divergence
    rho_conv  = drhou_x./(dx)  + drhov_y./(dy)  + drhow_z./(dz);
    rhou_conv = drhouu_x./(dx) + drhouv_y./(dy) + drhouw_z./(dz);
    rhov_conv = drhouv_x./(dx) + drhovv_y./(dy) + drhovw_z./(dz);
    rhow_conv = drhouw_x./(dx) + drhovw_y./(dy) + drhoww_z./(dz);
    rhoE_conv = drhouE_x./(dx) + drhovE_y./(dy) + drhowE_z./(dz);
    
    % Set outer points to 0
    rho_conv([1,end],:,:)    = 0;
    rho_conv(:,[1,end],:)    = 0;
    rho_conv(:,:,[1,end])    = 0;
    rhou_conv([1,end],:,:)   = 0;
    rhou_conv(:,[1,end],:)   = 0;
    rhou_conv(:,:,[1,end])   = 0;
    rhov_conv([1,end],:,:)   = 0;
    rhov_conv(:,[1,end],:)   = 0;
    rhov_conv(:,:,[1,end])   = 0;
    rhow_conv([1,end],:,:)   = 0;
    rhow_conv(:,[1,end],:)   = 0;
    rhow_conv(:,:,[1,end])   = 0;
    rhoE_conv([1,end],:,:)   = 0;
    rhoE_conv(:,[1,end],:)   = 0;
    rhoE_conv(:,:,[1,end])   = 0;

    C_D = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv};

    % C_u - Advective
    rho_conv  = u.*drho_x./(dx)  + v.*drho_y./(dy)  + w.*drho_z./(dz) + ...
         rho.*du_x./(dx)  + rho.*dv_y./(dy)  + rho.*dw_z./(dz);
    rhou_conv = u.*drhou_x./(dx)  + v.*drhou_y./(dy)  + w.*drhou_z./(dz) + ...
         rhou.*du_x./(dx)  + rhou.*dv_y./(dy)  + rhou.*dw_z./(dz);
    rhov_conv = u.*drhov_x./(dx)  + v.*drhov_y./(dy)  + w.*drhov_z./(dz) + ...
        rhov.*du_x./(dx)  + rhov.*dv_y./(dy)  + rhov.*dw_z./(dz);
    rhow_conv = u.*drhow_x./(dx)  + v.*drhow_y./(dy)  + w.*drhow_z./(dz) + ...
        rhow.*du_x./(dx)  + rhow.*dv_y./(dy)  + rhow.*dw_z./(dz);
    rhoE_conv =  u.*drhoE_x./(dx)  + v.*drhoE_y./(dy)  + w.*drhoE_z./(dz) + ...
        rhoE.*du_x./(dx)  + rhoE.*dv_y./(dy)  + rhoE.*dw_z./(dz);
        
    % Set outer points to 0
    rho_conv([1,end],:,:)    = 0;
    rho_conv(:,[1,end],:)    = 0;
    rho_conv(:,:,[1,end])    = 0;
    rhou_conv([1,end],:,:)   = 0;
    rhou_conv(:,[1,end],:)   = 0;
    rhou_conv(:,:,[1,end])   = 0;
    rhov_conv([1,end],:,:)   = 0;
    rhov_conv(:,[1,end],:)   = 0;
    rhov_conv(:,:,[1,end])   = 0;
    rhow_conv([1,end],:,:)   = 0;
    rhow_conv(:,[1,end],:)   = 0;
    rhow_conv(:,:,[1,end])   = 0;
    rhoE_conv([1,end],:,:)   = 0;
    rhoE_conv(:,[1,end],:)   = 0;
    rhoE_conv(:,:,[1,end])   = 0;
    
    C_u = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv};

    % C_fi - Advective
    rho_conv  = drhou_x./(dx)  + drhov_y./(dy)  + drhow_z./(dz) + ...
        0*rhou./(dx)  + 0*rhov./(dy)  + 0*rhow./(dz);
    rhou_conv = u.*drhou_x./(dx)  + u.*drhov_y./(dy)  + u.*drhow_z./(dz) + ...
        rhou.*du_x./(dx)  + rhov.*du_y./(dy)  + rhow.*du_z./(dz);
    rhov_conv = v.*drhou_x./(dx)  + v.*drhov_y./(dy)  + v.*drhow_z./(dz) + ...
        rhou.*dv_x./(dx)  + rhov.*dv_y./(dy)  + rhow.*dv_z./(dz);
    rhow_conv = w.*drhou_x./(dx)  + w.*drhov_y./(dy)  + w.*drhow_z./(dz) + ...
        rhou.*dw_x./(dx)  + rhov.*dw_y./(dy)  + rhow.*dw_z./(dz);
    rhoE_conv =  E.*drhou_x./(dx)  + E.*drhov_y./(dy)  + E.*drhow_z./(dz) + ...
        rhou.*dE_x./(dx)  + rhov.*dE_y./(dy)  + rhow.*dE_z./(dz);

    % Set outer points to 0
    rho_conv([1,end],:,:)    = 0;
    rho_conv(:,[1,end],:)    = 0;
    rho_conv(:,:,[1,end])    = 0;
    rhou_conv([1,end],:,:)   = 0;
    rhou_conv(:,[1,end],:)   = 0;
    rhou_conv(:,:,[1,end])   = 0;
    rhov_conv([1,end],:,:)   = 0;
    rhov_conv(:,[1,end],:)   = 0;
    rhov_conv(:,:,[1,end])   = 0;
    rhow_conv([1,end],:,:)   = 0;
    rhow_conv(:,[1,end],:)   = 0;
    rhow_conv(:,:,[1,end])   = 0;
    rhoE_conv([1,end],:,:)   = 0;
    rhoE_conv(:,[1,end],:)   = 0;
    rhoE_conv(:,:,[1,end])   = 0;
    
    C_fi = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv};


    % C_rho - Advective KG
    rho_conv  = rho.*du_x./(dx)  + rho.*dv_y./(dy)  + rho.*dw_z./(dz) + ...
        u.*drho_x./(dx)  + v.*drho_y./(dy)  + w.*drho_z./(dz);
    rhou_conv = rho.*duu_x./(dx)  + rho.*duv_y./(dy)  + rho.*duw_z./(dz) + ...
        u.*u.*drho_x./(dx)  + u.*v.*drho_y./(dy)  + u.*w.*drho_z./(dz);
    rhov_conv = rho.*duv_x./(dx)  + rho.*dvv_y./(dy)  + rho.*dvw_z./(dz) + ...
        v.*u.*drho_x./(dx)  + v.*v.*drho_y./(dy)  + v.*w.*drho_z./(dz);
    rhow_conv = rho.*duw_x./(dx)  + rho.*dvw_y./(dy)  + rho.*dww_z./(dz) + ...
        w.*u.*drho_x./(dx)  + w.*v.*drho_y./(dy)  + w.*w.*drho_z./(dz);
    rhoE_conv =  rho.*duE_x./(dx)  + rho.*dvE_y./(dy)  + rho.*dwE_z./(dz) + ...
        E.*u.*drho_x./(dx)  + E.*v.*drho_y./(dy)  + E.*w.*drho_z./(dz);

    % Set outer points to 0
    rho_conv([1,end],:,:)    = 0;
    rho_conv(:,[1,end],:)    = 0;
    rho_conv(:,:,[1,end])    = 0;
    rhou_conv([1,end],:,:)   = 0;
    rhou_conv(:,[1,end],:)   = 0;
    rhou_conv(:,:,[1,end])   = 0;
    rhov_conv([1,end],:,:)   = 0;
    rhov_conv(:,[1,end],:)   = 0;
    rhov_conv(:,:,[1,end])   = 0;
    rhow_conv([1,end],:,:)   = 0;
    rhow_conv(:,[1,end],:)   = 0;
    rhow_conv(:,:,[1,end])   = 0;
    rhoE_conv([1,end],:,:)   = 0;
    rhoE_conv(:,[1,end],:)   = 0;
    rhoE_conv(:,:,[1,end])   = 0;
    
    C_rho = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv};


    %% C_L - Advective linear
    rho_conv  = rho.*du_x./(dx)  + rho.*dv_y./(dy)  + rho.*dw_z./(dz) + ...
        0*rho.*u./(dx)  + 0*rho.*v./(dy)  + 0*rho.*w./(dz) + ...
        u.*drho_x./(dx)    + v.*drho_y./(dy)      + w.*drho_z./(dz);
    rhou_conv  = rho.*u.*du_x./(dx)  + rho.*u.*dv_y./(dy)  + rho.*u.*dw_z./(dz) + ...
        rho.*u.*du_x./(dx)  + rho.*v.*du_y./(dy)  + rho.*w.*du_z./(dz) + ...
        u.*u.*drho_x./(dx)  + u.*v.*drho_y./(dy)  + u.*w.*drho_z./(dz);
    rhov_conv  = rho.*v.*du_x./(dx)  + rho.*v.*dv_y./(dy)  + rho.*v.*dw_z./(dz) + ...
        rho.*u.*dv_x./(dx)  + rho.*v.*dv_y./(dy)  + rho.*w.*dv_z./(dz) + ...
        v.*u.*drho_x./(dx)  + v.*v.*drho_y./(dy)  + v.*w.*drho_z./(dz);
    rhow_conv  = rho.*w.*du_x./(dx)  + rho.*w.*dv_y./(dy)  + rho.*w.*dw_z./(dz) + ...
        rho.*u.*dw_x./(dx)  + rho.*v.*dw_y./(dy)  + rho.*w.*dw_z./(dz) + ...
        w.*u.*drho_x./(dx)  + w.*v.*drho_y./(dy)  + w.*w.*drho_z./(dz);
    rhoE_conv  = rho.*E.*du_x./(dx)  + rho.*E.*dv_y./(dy)  + rho.*E.*dw_z./(dz) + ...
        rho.*u.*dE_x./(dx)  + rho.*v.*dE_y./(dy)  + rho.*w.*dE_z./(dz) + ...
        E.*u.*drho_x./(dx)  + E.*v.*drho_y./(dy)  + E.*w.*drho_z./(dz);

    % Set outer points to 0
    rho_conv([1,end],:,:)    = 0;
    rho_conv(:,[1,end],:)    = 0;
    rho_conv(:,:,[1,end])    = 0;
    rhou_conv([1,end],:,:)   = 0;
    rhou_conv(:,[1,end],:)   = 0;
    rhou_conv(:,:,[1,end])   = 0;
    rhov_conv([1,end],:,:)   = 0;
    rhov_conv(:,[1,end],:)   = 0;
    rhov_conv(:,:,[1,end])   = 0;
    rhow_conv([1,end],:,:)   = 0;
    rhow_conv(:,[1,end],:)   = 0;
    rhow_conv(:,:,[1,end])   = 0;
    rhoE_conv([1,end],:,:)   = 0;
    rhoE_conv(:,[1,end],:)   = 0;
    rhoE_conv(:,:,[1,end])   = 0;

    C_L = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv};

end

