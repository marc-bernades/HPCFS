function [rho, rhou, rhov, rhow, rhoE] = Update_Boundaries_Simple(rho,rhou, rhov, rhow, rhoE, BC)

switch BC
    case 0 % Periodic
        % West boundary points
        rho(1,:,:)     = rho(end-1,:,:);
        rhou(1,:,:)    = rhou(end-1,:,:);
        rhov(1,:,:)    = rhov(end-1,:,:);
        rhow(1,:,:)    = rhow(end-1,:,:);
        rhoE(1,:,:)    = rhoE(end-1,:,:);

        % East boundary points
        rho(end,:,:)   = rho(2,:,:);
        rhou(end,:,:)  = rhou(2,:,:);
        rhov(end,:,:)  = rhov(2,:,:);
        rhow(end,:,:)  = rhow(2,:,:);
        rhoE(end,:,:)  = rhoE(2,:,:);

        % South boundary points
        rho(:,1,:)     = rho(:,end-1,:);
        rhou(:,1,:)    = rhou(:,end-1,:);
        rhov(:,1,:)    = rhov(:,end-1,:);
        rhow(:,1,:)    = rhow(:,end-1,:);
        rhoE(:,1,:)    = rhoE(:,end-1,:);

        % North boundary points
        rho(:,end,:)   = rho(:,2,:);
        rhou(:,end,:)  = rhou(:,2,:);
        rhov(:,end,:)  = rhov(:,2,:);
        rhow(:,end,:)  = rhow(:,2,:);
        rhoE(:,end,:)  = rhoE(:,2,:);

        % Back boundary points
        rho(:,:,1)     = rho(:,:,end-1);
        rhou(:,:,1)    = rhou(:,:,end-1);
        rhov(:,:,1)    = rhov(:,:,end-1);
        rhow(:,:,1)    = rhow(:,:,end-1);
        rhoE(:,:,1)    = rhoE(:,:,end-1);

        % Front boundary points
        rho(:,:,end)   = rho(:,:,2);
        rhou(:,:,end)  = rhou(:,:,2);
        rhov(:,:,end)  = rhov(:,:,2);
        rhow(:,:,end)  = rhow(:,:,2);
        rhoE(:,:,end)  = rhoE(:,:,2);

    case 1 % Neumann
        % West boundary points
        rho(1,:,:)     = rho(2,:,:);
        rhou(1,:,:)    = rhou(2,:,:);
        rhov(1,:,:)    = rhov(2,:,:);
        rhow(1,:,:)    = rhow(2,:,:);
        rhoE(1,:,:)    = rhoE(2,:,:);

        % East boundary points
        rho(end,:,:)   = rho(end-1,:,:);
        rhou(end,:,:)  = rhou(end-1,:,:);
        rhov(end,:,:)  = rhov(end-1,:,:);
        rhow(end,:,:)  = rhow(end-1,:,:);
        rhoE(end,:,:)  = rhoE(end-1,:,:);

        % South boundary points
        rho(:,1,:)     = rho(:,2,:);
        rhou(:,1,:)    = rhou(:,2,:);
        rhov(:,1,:)    = rhov(:,2,:);
        rhow(:,1,:)    = rhow(:,2,:);
        rhoE(:,1,:)    = rhoE(:,2,:);

        % North boundary points
        rho(:,end,:)   = rho(:,end-1,:);
        rhou(:,end,:)  = rhou(:,end-1,:);
        rhov(:,end,:)  = rhov(:,end-1,:);
        rhow(:,end,:)  = rhow(:,end-1,:);
        rhoE(:,end,:)  = rhoE(:,end-1,:);

        % Back boundary points
        rho(:,:,1)     = rho(:,:,2);
        rhou(:,:,1)    = rhou(:,:,2);
        rhov(:,:,1)    = rhov(:,:,2);
        rhow(:,:,1)    = rhow(:,:,2);
        rhoE(:,:,1)    = rhoE(:,:,2);

        % Front boundary points
        rho(:,:,end)   = rho(:,:,end-1);
        rhou(:,:,end)  = rhou(:,:,end-1);
        rhov(:,:,end)  = rhov(:,:,end-1);
        rhow(:,:,end)  = rhow(:,:,end-1);
        rhoE(:,:,end)  = rhoE(:,:,end-1);

    otherwise
end

 
end

