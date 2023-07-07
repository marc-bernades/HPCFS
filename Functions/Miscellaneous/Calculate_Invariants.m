function    [Invariants,Ref] = Calculate_Invariants(dx,dy,dz,...
        rho, u, v, w, E, e, c_v , gamma, T, P,t,Invariants,Ref,Substance,Test,bSolver)

% rho_norm_i, rhou_bar_i, rhoE_norm_i, rhoe_norm_i, rhos_norm_i, rhoui2_norm_i, ...
%         rhoE2_norm_i, rhoe2_norm_i, rhos2_norm_i,T_norm_i, e_norm_i, P_norm_i,

    % Volume
    volume        = dx.*dy.*dz;
    volume_sum    = sum(sum(sum(volume(2:end-1,2:end-1,2:end-1))));

    
    % Entropy calculations - If it is 1D test don't do anything on entropy
    % as it is imaginary
    if strcmp(Test,'1D_Adv_Shima') && strcmp(bSolver,'Real')
        rhos    = 1.0 + zeros(size(volume_sum));
        s      = rhos;
        sPR    = rhos;
        rhosPR = rhos;

    else
        % Calculate rho·s > Ideal gas relationship (ref JCP 2019 Coppola et
        % al.)
        rhos = rho.*c_v.*log((gamma - 1).*(rho.*E - 0.5*rho.*(u.^2 + v.^2 + w.^2))./(rho.^(gamma - 1)));
        s = rhos./rho;

        % Calculate s > Real gas framework
        sPR    = Calculate_Entropy(T,rho,P,Substance);
        rhosPR = rho.*sPR;

    end

    % Volumetric quantities
    rho_volume    = rho.*volume;
    rhou_volume   = rho.*u.*volume;
    rhoE_volume   = rho.*E.*volume;
    rhoe_volume   = rho.*e.*volume;
    rhos_volume   = rhos.*volume;
    rhosPR_volume = rhosPR.*volume;
    s_volume      = s.*volume;
    sPR_volume    = sPR.*volume;
    rhoui2_volume = rho.*(u.^2 + v.^2 + w.^2).*volume;
    rhoE2_volume  = rho.*E.^2.*volume;
    rhoe2_volume  = rho.*e.^2.*volume;
    rhos2_volume  = rhos.^2./rho.*volume;
    rhosPR2_volume= rhosPR.^2./rho.*volume;
    T_volume      = T.*volume;
    e_volume      = e.*volume;
    P_volume      = P.*volume;


    % Sum across the domain
    rho_volume_sum    = sum(sum(sum(rho_volume(2:end-1,2:end-1,2:end-1))));
    rhou_volume_sum   = sum(sum(sum(rhou_volume(2:end-1,2:end-1,2:end-1))));
    rhoE_volume_sum   = sum(sum(sum(rhoE_volume(2:end-1,2:end-1,2:end-1))));
    rhoe_volume_sum   = sum(sum(sum(rhoe_volume(2:end-1,2:end-1,2:end-1))));
    rhos_volume_sum   = sum(sum(sum(rhos_volume(2:end-1,2:end-1,2:end-1))));
    rhosPR_volume_sum = sum(sum(sum(rhosPR_volume(2:end-1,2:end-1,2:end-1))));
    s_volume_sum      = sum(sum(sum(s_volume(2:end-1,2:end-1,2:end-1))));
    sPR_volume_sum    = sum(sum(sum(sPR_volume(2:end-1,2:end-1,2:end-1))));
    rhoui2_volume_sum = sum(sum(sum(rhoui2_volume(2:end-1,2:end-1,2:end-1))));
    rhoE2_volume_sum  = sum(sum(sum(rhoE2_volume(2:end-1,2:end-1,2:end-1))));
    rhoe2_volume_sum  = sum(sum(sum(rhoe2_volume(2:end-1,2:end-1,2:end-1))));
    rhos2_volume_sum  = sum(sum(sum(rhos2_volume(2:end-1,2:end-1,2:end-1))));
    rhosPR2_volume_sum= sum(sum(sum(rhosPR2_volume(2:end-1,2:end-1,2:end-1))));
    T_volume_sum      = sum(sum(sum(T_volume(2:end-1,2:end-1,2:end-1))));
    e_volume_sum      = sum(sum(sum(e_volume(2:end-1,2:end-1,2:end-1))));
    P_volume_sum      = sum(sum(sum(P_volume(2:end-1,2:end-1,2:end-1))));

    % Total
    rho_norm_i     = rho_volume_sum/volume_sum;
    rhou_bar_i     = rhou_volume_sum/volume_sum;
    rhoE_norm_i    = rhoE_volume_sum/volume_sum;
    rhoe_norm_i    = rhoe_volume_sum/volume_sum;
    rhos_norm_i    = rhos_volume_sum/volume_sum;
    rhosPR_norm_i  = rhosPR_volume_sum/volume_sum;
    s_norm_i       = s_volume_sum/volume_sum;
    sPR_norm_i     = sPR_volume_sum/volume_sum;
    rhoui2_norm_i  = rhoui2_volume_sum/volume_sum;
    rhoE2_norm_i   = rhoE2_volume_sum/volume_sum;
    rhoe2_norm_i   = rhoe2_volume_sum/volume_sum;
    rhos2_norm_i   = rhos2_volume_sum/volume_sum;
    rhosPR2_norm_i = rhosPR2_volume_sum/volume_sum;
    T_norm_i       = T_volume_sum/volume_sum;    
    e_norm_i       = e_volume_sum/volume_sum;
    P_norm_i       = P_volume_sum/volume_sum;


    if t == 1
        Ref.rho_norm_ref    = rho_norm_i;
        Ref.rhoE_norm_ref   = rhoE_norm_i;
        Ref.rhoe_norm_ref   = rhoe_norm_i;
        Ref.rhos_norm_ref   = rhos_norm_i;
        Ref.rhosPR_norm_ref = rhosPR_norm_i;
        Ref.s_norm_ref      = s_norm_i;
        Ref.sPR_norm_ref    = sPR_norm_i;
        Ref.rhoui2_norm_ref = rhoui2_norm_i;
        Ref.rhoE2_norm_ref  = rhoE2_norm_i;
        Ref.rhoe2_norm_ref  = rhoe2_norm_i;
        Ref.rhos2_norm_ref  = rhos2_norm_i;
        Ref.rhosPR2_norm_ref= rhosPR2_norm_i;
        Ref.T_norm_ref      = T_norm_i;
        Ref.e_norm_ref      = e_norm_i;
        Ref.P_norm_ref      = P_norm_i;
    end
    Invariants.rho_norm(t)      = (rho_norm_i - Ref.rho_norm_ref)/abs(Ref.rho_norm_ref);
    Invariants.rhou_bar(t)      = rhou_bar_i;
    Invariants.rhoE_norm(t)     = (rhoE_norm_i - Ref.rhoE_norm_ref)/abs(Ref.rhoE_norm_ref);
    Invariants.rhoe_norm(t)     = (rhoe_norm_i - Ref.rhoe_norm_ref)/abs(Ref.rhoe_norm_ref);
    Invariants.rhos_norm(t)     = (rhos_norm_i - Ref.rhos_norm_ref)/abs(Ref.rhos_norm_ref);
    Invariants.rhosPR_norm(t)   = (rhosPR_norm_i - Ref.rhosPR_norm_ref)/abs(Ref.rhosPR_norm_ref);
    Invariants.s_norm(t)        = (s_norm_i - Ref.s_norm_ref)/abs(Ref.s_norm_ref);
    Invariants.sPR_norm(t)      = (sPR_norm_i - Ref.sPR_norm_ref)/abs(Ref.sPR_norm_ref);
    Invariants.rhoui2_norm(t)   = (rhoui2_norm_i - Ref.rhoui2_norm_ref)/abs(Ref.rhoui2_norm_ref);
    Invariants.rhoE2_norm(t)    = (rhoE2_norm_i - Ref.rhoE2_norm_ref)/abs(Ref.rhoE2_norm_ref);
    Invariants.rhoe2_norm(t)    = (rhoe2_norm_i - Ref.rhoe2_norm_ref)/abs(Ref.rhoe2_norm_ref);
    Invariants.rhos2_norm(t)    = (rhos2_norm_i - Ref.rhos2_norm_ref)/abs(Ref.rhos2_norm_ref);
    Invariants.rhosPR2_norm(t)  = (rhosPR2_norm_i - Ref.rhosPR2_norm_ref)/abs(Ref.rhosPR2_norm_ref);
    Invariants.T_norm(t)        = (T_norm_i - Ref.T_norm_ref)/abs(Ref.T_norm_ref);
    Invariants.e_norm(t)        = (e_norm_i - Ref.e_norm_ref)/abs(Ref.e_norm_ref);
    Invariants.P_norm(t)        = (P_norm_i - Ref.P_norm_ref)/abs(Ref.P_norm_ref);



end