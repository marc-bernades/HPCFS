function [C_D, C_u, C_fi, C_rho, C_L] = ConvectiveTermSplit(rho,u,v,w,E,X,Y,Z,P)
    
    % Delta x, y and z based on CentralDerivative_d1_2ndOrder
    [dx,~,~] = CentralDerivative_d1_2ndOrder(X);
    [~,dy,~] = CentralDerivative_d1_2ndOrder(Y);
    [~,~,dz] = CentralDerivative_d1_2ndOrder(Z);

    % Calculate enthalpy to make its splitting (includes pressure CentralDerivative_d1_2ndOrder)
    k    = 0.5.*(u.^2 + v.^2 + w.^2);
    e    = E - k; % Comment for rhoe evolution
%     h    = E + P./rho;

    %% Convective terms
    % CD - Divergence
    rho_conv  = Inviscid_CD(rho,u,v,w,ones(size(rho)),dx,dy,dz);
    rhou_conv = Inviscid_CD(rho,u,v,w,u,dx,dy,dz);
    rhov_conv = Inviscid_CD(rho,u,v,w,v,dx,dy,dz);
    rhow_conv = Inviscid_CD(rho,u,v,w,w,dx,dy,dz);
    rhoe_conv = Inviscid_CD(rho,u,v,w,e,dx,dy,dz);
    rhok_conv = Inviscid_CD(rho,u,v,w,k,dx,dy,dz);
    Pu_conv   = Inviscid_CD(P,u,v,w,ones(size(rho)),dx,dy,dz);

    C_D = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoe_conv, rhok_conv, Pu_conv};

    % C_u - Advective
    rho_conv  = Inviscid_Cu(rho,u,v,w,ones(size(rho)),dx,dy,dz);
    rhou_conv = Inviscid_Cu(rho,u,v,w,u,dx,dy,dz);
    rhov_conv = Inviscid_Cu(rho,u,v,w,v,dx,dy,dz);
    rhow_conv = Inviscid_Cu(rho,u,v,w,w,dx,dy,dz);
    rhoe_conv = Inviscid_Cu(rho,u,v,w,e,dx,dy,dz);
    rhok_conv = Inviscid_Cu(rho,u,v,w,k,dx,dy,dz);
    Pu_conv   = Inviscid_Cu(P,u,v,w,ones(size(rho)),dx,dy,dz);

    C_u = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoe_conv, rhok_conv, Pu_conv};

    % C_fi - Advective
    rho_conv  = Inviscid_Cfi(rho,u,v,w,ones(size(rho)),dx,dy,dz);
    rhou_conv = Inviscid_Cfi(rho,u,v,w,u,dx,dy,dz);
    rhov_conv = Inviscid_Cfi(rho,u,v,w,v,dx,dy,dz);
    rhow_conv = Inviscid_Cfi(rho,u,v,w,w,dx,dy,dz);
    rhoe_conv = Inviscid_Cfi(rho,u,v,w,e,dx,dy,dz);
    rhok_conv = Inviscid_Cfi(rho,u,v,w,k,dx,dy,dz);
    Pu_conv   = Inviscid_Cfi(P,u,v,w,ones(size(rho)),dx,dy,dz);

    C_fi = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoe_conv, rhok_conv, Pu_conv};

    % C_rho - Advective KG
    rho_conv  = Inviscid_Crho(rho,u,v,w,ones(size(rho)),dx,dy,dz);
    rhou_conv = Inviscid_Crho(rho,u,v,w,u,dx,dy,dz);
    rhov_conv = Inviscid_Crho(rho,u,v,w,v,dx,dy,dz);
    rhow_conv = Inviscid_Crho(rho,u,v,w,w,dx,dy,dz);
    rhoe_conv = Inviscid_Crho(rho,u,v,w,e,dx,dy,dz);
    rhok_conv = Inviscid_Crho(rho,u,v,w,k,dx,dy,dz);
    Pu_conv   = Inviscid_Crho(P,u,v,w,ones(size(rho)),dx,dy,dz);

    C_rho = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoe_conv, rhok_conv, Pu_conv};

    %% C_L - Advective linear
    rho_conv  = Inviscid_CL(rho,u,v,w,ones(size(rho)),dx,dy,dz);
    rhou_conv = Inviscid_CL(rho,u,v,w,u,dx,dy,dz);
    rhov_conv = Inviscid_CL(rho,u,v,w,v,dx,dy,dz);
    rhow_conv = Inviscid_CL(rho,u,v,w,w,dx,dy,dz);
    rhoe_conv = Inviscid_CL(rho,u,v,w,e,dx,dy,dz);
    rhok_conv = Inviscid_CL(rho,u,v,w,k,dx,dy,dz);
    Pu_conv   = Inviscid_CL(P,u,v,w,ones(size(rho)),dx,dy,dz);

    C_L = {rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoe_conv, rhok_conv, Pu_conv}; 

end

