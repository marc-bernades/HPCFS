function [rho_inv_flux,rhou_inv_flux,rhov_inv_flux,rhow_inv_flux,rhoE_inv_flux] = CalculateInviscidFlux_WENO5_Projection(rho,u,v,w,E,sos,P,dx,dy,dz,bFlux,bSolver, Fluid, Substance)

% Initialize
rho_inv_flux  = zeros(size(u));
rhou_inv_flux = zeros(size(u));
rhov_inv_flux = zeros(size(u));
rhow_inv_flux = zeros(size(u));
rhoE_inv_flux = zeros(size(u));

Flux_x        = zeros([size(u),5]);
Flux_y        = zeros([size(u),5]);
Flux_z        = zeros([size(u),5]);

% Compute conserved
rho    = rho;
rhou   = rho.*u;
rhov   = rho.*v;
rhow   = rho.*w;
rhoE   = rho.*E;

% Compute Roe Average flux at 1/2
for dimension = 1:3
    if dimension == 1
        u_dim = u;
        dx_input = dx;
        sym{1} = 1; sym{2} = 0; sym{3} = 0;
    elseif dimension == 2
        u_dim = v;
        dx_input = dy;
        sym{1} = 0; sym{2} = 1; sym{3} = 0;
    elseif dimension == 3
        u_dim = w;
        dx_input = dz;
        sym{1} = 0; sym{2} = 0; sym{3} = 1;
    end
    % Conserved
    C{dimension}.rho    = Calculate_Roe_Avg(rho,dx_input,dimension);
    C{dimension}.rhou   = Calculate_Roe_Avg(rhou,dx_input,dimension);
    C{dimension}.rhov   = Calculate_Roe_Avg(rhov,dx_input,dimension);
    C{dimension}.rhow   = Calculate_Roe_Avg(rhow,dx_input,dimension);
    C{dimension}.rhoE   = Calculate_Roe_Avg(rhoE,dx_input,dimension);
    C{dimension}.u_dim  = Calculate_Roe_Avg(u_dim,dx_input,dimension);
    C{dimension}.P      = Calculate_Roe_Avg(P,dx_input,dimension);
    % Flux
    F{dimension}.rho    = C{dimension}.rho;
    F{dimension}.rhou   = C{dimension}.rho + sym{1}*C{dimension}.P;
    F{dimension}.rhov   = C{dimension}.rho + sym{2}*C{dimension}.P;
    F{dimension}.rhow   = C{dimension}.rho + sym{3}*C{dimension}.P;
    F{dimension}.rhoE   = C{dimension}.rho + C{dimension}.P.*C{dimension}.u_input;
end

% Projection to characteristic space
% Point-wise calculation
for i = 2:length(u(:,1,1))-1
    for j = 2:length(u(1,:,1))-1
        for k = 2:length(u(1,1,:))-1
            for d = 1:3
                % Calculate L matrix
                L = Calculate_ProjectionMatrix(F,C,dimension,i,j,k)

            end

        end
    end
end

% Only inner points
for i = 2:length(u(:,1,1))-1
    for j = 2:length(u(1,:,1))-1
        for k = 2:length(u(1,1,:))-1
            
            %% X - Direction (j component)
            dimension = 1;
            if length(u(1,:,1)) > 8 % Stencil 6 points + 2 boundaries = 8
                % Projection to charateristic space
                
                [Fp,Fn] = WENO5_Stencil_Calculation(i,j,k,rho,u,v,w,E,P,sos, dimension);
                % Projection to characteristic space

                % WENO 5 > Flux for each of the variables 1:5
                Flux_x(i,j,k,:)  = WENO5_solver(Fp,Fn);
            end

            %% Y - Direction (i component)
            dimension = 2;
            if length(u(:,1,1)) > 8 % Stencil 6 points + 2 boundaries = 8
                % Calculate fl ux i +/- 1/2
                [Fp,Fn] = WENO5_Stencil_Calculation(i,j,k,rho,u,v,w,E,P,sos, dimension);
                % WENO 5 > Flux for each of the variables 1:5
                Flux_y(i,j,k,:)  = WENO5_solver(Fp,Fn);
            end

            %% Z - Direction (k component)        
            dimension = 3;
            if length(u(1,1,:)) > 8 % Stencil 6 points + 2 boundaries = 8
                % Calculate flux i +/- 1/2
                [Fp,Fn] = WENO5_Stencil_Calculation(i,j,k,rho,u,v,w,E,P,sos, dimension);
                % WENO 5 > Flux for each of the variables 1:5
                Flux_z(i,j,k,:)  = WENO5_solver(Fp,Fn);
            end

        end
    end
end

% Update first point flux (last not needed) to have flux at 1/2
Flux_x(:,1,:,1)   = Flux_x(:,end-1,:,1);
Flux_x(:,1,:,2)   = Flux_x(:,end-1,:,2);
Flux_x(:,1,:,3)   = Flux_x(:,end-1,:,3);
Flux_x(:,1,:,4)   = Flux_x(:,end-1,:,4);
Flux_x(:,1,:,5)   = Flux_x(:,end-1,:,5);

Flux_y(1,:,:,1)   = Flux_y(end-1,:,:,1);
Flux_y(1,:,:,2)   = Flux_y(end-1,:,:,2);
Flux_y(1,:,:,3)   = Flux_y(end-1,:,:,3);
Flux_y(1,:,:,4)   = Flux_y(end-1,:,:,4);
Flux_y(1,:,:,5)   = Flux_y(end-1,:,:,5);

Flux_z(:,:,1,1)   = Flux_z(:,:,end-1,1);
Flux_z(:,:,1,2)   = Flux_z(:,:,end-1,2);
Flux_z(:,:,1,3)   = Flux_z(:,:,end-1,3);
Flux_z(:,:,1,4)   = Flux_z(:,:,end-1,4);
Flux_z(:,:,1,5)   = Flux_z(:,:,end-1,5);


% Only inner points
for i = 2:length(u(:,1,1))-1
    for j = 2:length(u(1,:,1))-1
        for k = 2:length(u(1,1,:))-1

            % Compute flux x-direction > x-dir (f_{i+1/2} - f_{i-1/2}) /dx
            rho_inv_flux(i,j,k)   = (Flux_x(i,j,k,1) - Flux_x(i,j-1,k,1))./dx(i,j,k);
            rhou_inv_flux(i,j,k)  = (Flux_x(i,j,k,2) - Flux_x(i,j-1,k,2))./dx(i,j,k);
            rhov_inv_flux(i,j,k)  = (Flux_x(i,j,k,3) - Flux_x(i,j-1,k,3))./dx(i,j,k);
            rhow_inv_flux(i,j,k)  = (Flux_x(i,j,k,4) - Flux_x(i,j-1,k,4))./dx(i,j,k);
            rhoE_inv_flux(i,j,k)  = (Flux_x(i,j,k,5) - Flux_x(i,j-1,k,5))./dx(i,j,k);

            % Compute flux y-direction
            rho_inv_flux(i,j,k)  = rho_inv_flux(i,j,k)  + (Flux_y(i,j,k,1) - Flux_y(i-1,j,k,1))./dy(i,j,k);
            rhou_inv_flux(i,j,k) = rhou_inv_flux(i,j,k) + (Flux_y(i,j,k,2) - Flux_y(i-1,j,k,2))./dy(i,j,k);
            rhov_inv_flux(i,j,k) = rhov_inv_flux(i,j,k) + (Flux_y(i,j,k,3) - Flux_y(i-1,j,k,3))./dy(i,j,k);
            rhow_inv_flux(i,j,k) = rhow_inv_flux(i,j,k) + (Flux_y(i,j,k,4) - Flux_y(i-1,j,k,4))./dy(i,j,k);
            rhoE_inv_flux(i,j,k) = rhoE_inv_flux(i,j,k) + (Flux_y(i,j,k,5) - Flux_y(i-1,j,k,5))./dy(i,j,k);

            % Compute flux z-direction
            rho_inv_flux(i,j,k)  = rho_inv_flux(i,j,k)  + (Flux_z(i,j,k,1) - Flux_z(i,j,k-1,1))./dz(i,j,k);
            rhou_inv_flux(i,j,k) = rhou_inv_flux(i,j,k) + (Flux_z(i,j,k,2) - Flux_z(i,j,k-1,2))./dz(i,j,k);
            rhov_inv_flux(i,j,k) = rhov_inv_flux(i,j,k) + (Flux_z(i,j,k,3) - Flux_z(i,j,k-1,3))./dz(i,j,k);
            rhow_inv_flux(i,j,k) = rhow_inv_flux(i,j,k) + (Flux_z(i,j,k,4) - Flux_z(i,j,k-1,4))./dz(i,j,k);
            rhoE_inv_flux(i,j,k) = rhoE_inv_flux(i,j,k) + (Flux_z(i,j,k,5) - Flux_z(i,j,k-1,5))./dz(i,j,k);


        end
    end
end


end