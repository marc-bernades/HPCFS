function [C_D, C_u, C_fi, C_rho, C_L] = ConvectiveTermShima(rho,u,v,w,E,rhou,rhov,rhow,rhoE,X,Y,Z,P)
    
    % Delta x, y and z based on CentralDerivative_d1_2ndOrder
    [dx,~,~] = CentralDerivative_d1_2ndOrder(X);
    [~,dy,~] = CentralDerivative_d1_2ndOrder(Y);
    [~,~,dz] = CentralDerivative_d1_2ndOrder(Z);

    % Total split in rho u h = rho u (e + ke + p/rho)
    h    = E + P./rho;
    rhoh = rho.*h;
    ke   = 0.5*(u.^2 + v.^2 + w.^2);
    e    = E - ke;

    % CentralDerivative_d1_2ndOrder computes de 2nd order derivative with respect to x,y,z
    [du_x,du_y,du_z]             = CentralDerivative_d1_2ndOrder(u);
    [dv_x,dv_y,dv_z]             = CentralDerivative_d1_2ndOrder(v);
    [dw_x,dw_y,dw_z]             = CentralDerivative_d1_2ndOrder(w);
    [drho_x,drho_y,drho_z]       = CentralDerivative_d1_2ndOrder(rho);
    [dh_x,  dh_y,  dh_z]         = CentralDerivative_d1_2ndOrder(h);

    [duu_x,~,~]                  = CentralDerivative_d1_2ndOrder(u.*u);
    [~,    dvv_y,~]              = CentralDerivative_d1_2ndOrder(v.*v);
    [~,    ~,    dww_z]          = CentralDerivative_d1_2ndOrder(w.*w);

    [duv_x,duv_y,~]              = CentralDerivative_d1_2ndOrder(u.*v);
    [duw_x,~,    duw_z]          = CentralDerivative_d1_2ndOrder(u.*w);
    [~,    dvw_y,dvw_z]          = CentralDerivative_d1_2ndOrder(v.*w);

    [duh_x,~,    ~]              = CentralDerivative_d1_2ndOrder(u.*h);
    [~,    dvh_y,~]              = CentralDerivative_d1_2ndOrder(v.*h);
    [~,    ~,    dwh_z]          = CentralDerivative_d1_2ndOrder(w.*h);

    [drhou_x,drhou_y,drhou_z]    = CentralDerivative_d1_2ndOrder(rhou);
    [drhov_x,drhov_y,drhov_z]    = CentralDerivative_d1_2ndOrder(rhov);
    [drhow_x,drhow_y,drhow_z]    = CentralDerivative_d1_2ndOrder(rhow);
    [drhoh_x,drhoh_y,drhoh_z]    = CentralDerivative_d1_2ndOrder(rhoh);

    [drhouu_x,drhouu_y,drhouu_z] = CentralDerivative_d1_2ndOrder(rhou.*u);
    [drhovv_x,drhovv_y,drhovv_z] = CentralDerivative_d1_2ndOrder(rhov.*v);
    [drhoww_x,drhoww_y,drhoww_z] = CentralDerivative_d1_2ndOrder(rhow.*w);

    [drhouv_x,drhouv_y,~]        = CentralDerivative_d1_2ndOrder(rhou.*v);
    [drhouw_x,~,       drhouw_z] = CentralDerivative_d1_2ndOrder(rhou.*w);
    [~,       drhovw_y,drhovw_z] = CentralDerivative_d1_2ndOrder(rhov.*w);
    
    [drhouh_x,~       ,~]        = CentralDerivative_d1_2ndOrder(rhou.*h);
    [~       ,drhovh_y,~]        = CentralDerivative_d1_2ndOrder(rhov.*h);
    [~       ,~       ,drhowh_z] = CentralDerivative_d1_2ndOrder(rhow.*h);

    % Prepare energy terms
    [drhoeu_x,drhoeu_y,drhoeu_z]  = CentralDerivative_d1_2ndOrder(rhou.*e.*u);
    [drhoev_x,drhoev_y,drhoev_z]  = CentralDerivative_d1_2ndOrder(rhou.*e.*v);
    [drhoew_x,drhoew_y,drhoew_z]  = CentralDerivative_d1_2ndOrder(rhou.*e.*w);

    [drhoe_x, drhoe_y, drhoe_z]   = CentralDerivative_d1_2ndOrder(rhou.*e);
    [dP_x,    dP_y,    dP_z]      = CentralDerivative_d1_2ndOrder(P);

    [de_x,de_y,de_z]  = CentralDerivative_d1_2ndOrder(e);

    [deu_x,deu_y,deu_z]  = CentralDerivative_d1_2ndOrder(e.*u);
    [dev_x,dev_y,dev_z]  = CentralDerivative_d1_2ndOrder(e.*v);
    [dew_x,dew_y,dew_z]  = CentralDerivative_d1_2ndOrder(e.*w);


    % Prepare ke terms
    [drhouuu_x,~,       ~]        = CentralDerivative_d1_2ndOrder(rhou.*u.*u);
    [~,       drhovvv_y,~]        = CentralDerivative_d1_2ndOrder(rhov.*v.*v);
    [~,        ~,      drhowww_z] = CentralDerivative_d1_2ndOrder(rhow.*w.*w);

    [drhouvv_x,drhouvv_y,drhouvv_z]  = CentralDerivative_d1_2ndOrder(rhou.*v.*v);
    [drhouww_x,drhouww_y,drhouww_z]  = CentralDerivative_d1_2ndOrder(rhou.*w.*w);
    [drhovuu_x,drhovuu_y,drhovuu_z]  = CentralDerivative_d1_2ndOrder(rhov.*u.*u);
    [drhovww_x,drhovww_y,drhovww_z]  = CentralDerivative_d1_2ndOrder(rhov.*w.*w);
    [drhowuu_x,drhowuu_y,drhowuu_z]  = CentralDerivative_d1_2ndOrder(rhow.*u.*u);
    [drhowvv_x,drhowvv_y,drhowvv_z]  = CentralDerivative_d1_2ndOrder(rhow.*v.*v);


%     rhouk_conv = 1/2*1/2*(drhouuu_x./dx + drhovuu_y./dy + drhowuu_z./dz + ...
%         drhouvv_x./dx + drhovvv_y./dy + drhowvv_z./dz + ...
%         drhouww_x./dx + drhovww_y./dy + drhowww_z./dz + ...
%         u.*(drhouu_x./dx + drhovv_x./dx + drhoww_x./dx) + ...
%         v.*(drhouu_y./dy + drhovv_y./dy + drhoww_y./dy) + ...
%         w.*(drhouu_z./dz + drhovv_z./dz + drhoww_z./dz) + ...
%         rho.*(u.^2 + v.^2 + w.^2).*du_x./dx + ... 
%         rho.*(u.^2 + v.^2 + w.^2).*dv_y./dy + ...
%         rho.*(u.^2 + v.^2 + w.^2).*dw_z./dz);

     rhouk_conv = 1/2*1/2*(u.*drhouu_x./dx + u.*drhouv_y./dy + u.*drhouw_z./dz + ...
        v.*drhouv_x./dx + v.*drhovv_y./dy + v.*drhovw_z./dz + ...
        w.*drhouw_x./dx + w.*drhovw_y./dy + w.*drhoww_z./dz + ...
        rho.*u.*u.*du_x./dx + rho.*u.*v.*du_y./dy + rho.*u.*w.*du_z./dz + ...
        rho.*v.*u.*dv_x./dx + rho.*v.*v.*dv_y./dy + rho.*v.*w.*dv_z./dz + ...
        rho.*w.*u.*dw_x./dx + rho.*w.*v.*dw_y./dy + rho.*w.*w.*dw_z./dz + ...
        u.*u.*drhou_x./dx + u.*v.*drhou_y./dy + u.*w.*drhou_z./dz + ...
        v.*u.*drhov_x./dx + v.*v.*drhov_y./dy + v.*w.*drhov_z./dz + ...
        w.*u.*drhow_x./dx + w.*v.*drhow_y./dy + w.*w.*drhow_z./dz + ...
        rho.*u.*drhouu_x./dx + rho.*u.*drhouv_y./dy + rho.*u.*drhouw_z./dz + ...
        rho.*v.*drhouv_x./dx + rho.*v.*drhovv_y./dy + rho.*v.*drhovw_z./dz + ...
        rho.*w.*drhouw_x./dx + rho.*w.*drhovw_y./dy + rho.*w.*drhoww_z./dz);


    %% Shima scheme
    % Internal energy
%     rhoeu_conv = 1/2*(drhoeu_x./dx + drhoeu_y./dy + drhoeu_z./dz + ...
%         u.*drhoe_x./dx + v.*drhoe_y./dy + w.*drhoe_z./dz + ...
%         rho.*e.*du_x./dx + rho.*e.*dv_y./dy + rho.*e.*dw_z./dz);

    rhoeu_conv = 1/4*(drhoeu_x./dx + drhoev_y./dy + drhoew_z./dz + ...
        u.*drhoe_x./dx + v.*drhoe_y./dy + w.*drhoe_z./dz + ...
        rho.*e.*(du_x./dx + dv_y./dy + dw_z./dz) + ...
        e.*(drhou_x./dx + drhov_y./dy + drhou_z./dz) + ...
        rho.*u.*de_x./dx + rho.*v.*de_y./dy + rho.*w.*de_z./dz + ...
        rho.*(deu_x./dx + dev_y./dy + dew_z./dz) + ...
        e.*u.*drho_x./dx +  e.*v.*drho_y./dy +  e.*w.*drho_z./dz);

    % Kinetic energy
%     rhouk_conv = 1/4*(drhouu_x./dx + drhovv_y./dy + drhoww_z./dz + ...
%         rho.*duu_x./dx + rho.*dvv_y./dy + rho.*dww_z./dz + ...
%         u.*drhouu_x./dx + v.*drhovv_y./dy + w.*drhoww_z./dz);
%     rhouk_conv_D = (drhouu_x./dx + drhouv_y./dy + drhouw_z./dz + ...
%         drhouv_x./dx + drhovv_y./dy + drhovw_z./dz + ...
%         drhouw_x./dx + drhovw_y./dy + drhoww_z./dz);
%     rhouk_conv_u = u.*drhou_x./(dx)  + v.*drhou_y./(dy)  + w.*drhou_z./(dz) + ...
%          rhou.*du_x./(dx)  + rhou.*dv_y./(dy)  + rhou.*dw_z./(dz) + ...
%          u.*drhov_x./(dx)  + v.*drhov_y./(dy)  + w.*drhov_z./(dz) + ...
%          rhov.*du_x./(dx)  + rhov.*dv_y./(dy)  + rhov.*dw_z./(dz) + ...
%          u.*drhow_x./(dx)  + v.*drhow_y./(dy)  + w.*drhow_z./(dz) + ...
%          rhow.*du_x./(dx)  + rhow.*dv_y./(dy)  + rhow.*dw_z./(dz);
%     rhouk_conv_fi = u.*drhou_x./(dx)  + u.*drhov_y./(dy)  + u.*drhow_z./(dz) + ...
%         rhou.*du_x./(dx)  + rhov.*du_y./(dy)  + rhow.*du_z./(dz) + ...
%         v.*drhou_x./(dx)  + v.*drhov_y./(dy)  + v.*drhow_z./(dz) + ...
%         rhou.*dv_x./(dx)  + rhov.*dv_y./(dy)  + rhow.*dv_z./(dz) + ...
%         w.*drhou_x./(dx)  + w.*drhov_y./(dy)  + w.*drhow_z./(dz) + ...
%         rhou.*dw_x./(dx)  + rhov.*dw_y./(dy)  + rhow.*dw_z./(dz);
%     rhouk_conv_rho = rho.*duu_x./(dx)  + rho.*duv_y./(dy)  + rho.*duw_z./(dz) + ...
%         u.*u.*drho_x./(dx)  + u.*v.*drho_y./(dy)  + u.*w.*drho_z./(dz) + ...
%         rho.*duv_x./(dx)  + rho.*dvv_y./(dy)  + rho.*dvw_z./(dz) + ...
%         v.*u.*drho_x./(dx)  + v.*v.*drho_y./(dy)  + v.*w.*drho_z./(dz) + ...
%         rho.*duw_x./(dx)  + rho.*dvw_y./(dy)  + rho.*dww_z./(dz) + ...
%         w.*u.*drho_x./(dx)  + w.*v.*drho_y./(dy)  + w.*w.*drho_z./(dz);
% 
%     rhouk_conv = 1/4*(rhouk_conv_D + rhouk_conv_u + rhouk_conv_fi + rhouk_conv_rho);
    
    % Pressure
%     rhouP_conv = 1/2*(u.*dP_x./dx + v.*dP_y./dy + w.*dP_z./dz + ...
%         P.*du_x./dx + P.*dv_y./dy + P.*dw_z./dz);

    rhouP_conv = u.*dP_x./dx + v.*dP_y./dy + w.*dP_z./dz + ...
        P.*du_x./dx + P.*dv_y./dy + P.*dw_z./dz;

    % Convective term for Shima scheme
    rhoE_conv  = rhoeu_conv + rhouk_conv + rhouP_conv;



    %% Convective terms
    % CD - Divergence
    rho_conv  = drhou_x./(dx)  + drhov_y./(dy)  + drhow_z./(dz);
    rhou_conv = drhouu_x./(dx) + drhouv_y./(dy) + drhouw_z./(dz);
    rhov_conv = drhouv_x./(dx) + drhovv_y./(dy) + drhovw_z./(dz);
    rhow_conv = drhouw_x./(dx) + drhovw_y./(dy) + drhoww_z./(dz);
%     rhoE_conv = drhouh_x./(dx) + drhovh_y./(dy) + drhowh_z./(dz);

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
%     rhoE_conv =  u.*drhoh_x./(dx)  + v.*drhoh_y./(dy)  + w.*drhoh_z./(dz) + ...
%         rhoh.*du_x./(dx)  + rhoh.*dv_y./(dy)  + rhoh.*dw_z./(dz);

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
%     rhoE_conv =  h.*drhou_x./(dx)  + h.*drhov_y./(dy)  + h.*drhow_z./(dz) + ...
%         rhou.*dh_x./(dx)  + rhov.*dh_y./(dy)  + rhow.*dh_z./(dz);

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
%     rhoE_conv =  rho.*duh_x./(dx)  + rho.*dvh_y./(dy)  + rho.*dwh_z./(dz) + ...
%         h.*u.*drho_x./(dx)  + h.*v.*drho_y./(dy)  + h.*w.*drho_z./(dz);

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
%     rhoE_conv  = rho.*h.*du_x./(dx)  + rho.*h.*dv_y./(dy)  + rho.*h.*dw_z./(dz) + ...
%         rho.*u.*dh_x./(dx)  + rho.*v.*dh_y./(dy)  + rho.*w.*dh_z./(dz) + ...
%         h.*u.*drho_x./(dx)  + h.*v.*drho_y./(dy)  + h.*w.*drho_z./(dz);

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

