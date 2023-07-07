function [rho, rhou, rhov, rhow, rhoE, P] = Update_Boundaries(bSolver,rho, rhou, rhov, rhow, rhoE, u, v, w, P, T, BC, Fluid,Substance)

Boundaries = fieldnames(BC);

for jj = 1:length(Boundaries)

    switch Boundaries{jj}
        case 'south'
            % Check type of boundary
            switch BC.(Boundaries{jj}){1}

                case 'Periodic'
                    % West boundary points
                    rho(1,:,:)     = rho(end-1,:,:);
                    rhou(1,:,:)    = rhou(end-1,:,:);
                    rhov(1,:,:)    = rhov(end-1,:,:);
                    rhow(1,:,:)    = rhow(end-1,:,:);
                    rhoE(1,:,:)    = rhoE(end-1,:,:);
                    P(1,:,:)       = P(end-1,:,:);


                case 'Neumann'
                    % West boundary points
                    rho(1,:,:)     = rho(2,:,:);
                    rhou(1,:,:)    = rhou(2,:,:);
                    rhov(1,:,:)    = -rhov(2,:,:);
                    rhow(1,:,:)    = rhow(2,:,:);
                    rhoE(1,:,:)    = rhoE(2,:,:);
                    P(1,:,:)       = P(2,:,:);

                case 'Dirichlet'
                    % West boundary points
                    u(1,:,:)       = 2*BC.(Boundaries{jj}){2} - u(2,:,:);
                    v(1,:,:)       = 2*BC.(Boundaries{jj}){3} - v(2,:,:);
                    w(1,:,:)       = 2*BC.(Boundaries{jj}){4} - w(2,:,:);
                    % Check impermeable boundary
                    if BC.(Boundaries{jj}){5} > 0
                        P(1,:,:)   = 2*BC.(Boundaries{jj}){5} - P(2,:,:);
                    elseif BC.(Boundaries{jj}){5} < 0
                        % Impermeable boundary
                        P(1,:,:)    = P(2,:,:);
                    else
                        % Don't do anything
                    end
                    % Check isothermal boundary
                    if BC.(Boundaries{jj}){6} > 0
                        T(1,:,:)        = 2*BC.(Boundaries{jj}){6} - T(2,:,:);
                    else
                        % No effect of temperature

                    end
                    % Update E and rho based on P and T
                    rho  = Calculate_Rho_from_TandP( bSolver, T,P, Fluid,Substance );
                    e    = Calculate_e_from_TandPandRho( bSolver, T,rho,P, Fluid, Substance);
                    ke  = 0.5*(u.^2 + v.^2 + w.^2);
                    E   = e + ke;
                    % Update conserved variables
                    rhou = rho.*u;
                    rhov = rho.*v;
                    rhow = rho.*w;
                    rhoE = rho.*E;

                otherwise

            end
        case 'north'
            % Check type of boundary
            switch BC.(Boundaries{jj}){1}

                case 'Periodic'
                    % East boundary points
                    rho(end,:,:)   = rho(2,:,:);
                    rhou(end,:,:)  = rhou(2,:,:);
                    rhov(end,:,:)  = rhov(2,:,:);
                    rhow(end,:,:)  = rhow(2,:,:);
                    rhoE(end,:,:)  = rhoE(2,:,:);
                    P(end,:,:)     = P(2,:,:);

                case 'Neumann'
                    % East boundary points
                    rho(end,:,:)   = rho(end-1,:,:);
                    rhou(end,:,:)  = rhou(end-1,:,:);
                    rhov(end,:,:)  = -rhov(end-1,:,:);
                    rhow(end,:,:)  = rhow(end-1,:,:);
                    rhoE(end,:,:)  = rhoE(end-1,:,:);
                    P(end,:,:)     = P(end-1,:,:);

                case 'Dirichlet'
                    % East boundary points
                    u(end,:,:)       = 2*BC.(Boundaries{jj}){2} - u(end-1,:,:);
                    v(end,:,:)       = 2*BC.(Boundaries{jj}){3} - v(end-1,:,:);
                    w(end,:,:)       = 2*BC.(Boundaries{jj}){4} - w(end-1,:,:);
                    % Check impermeable boundary
                    if BC.(Boundaries{jj}){5} > 0
                        P(end,:,:)   = 2*BC.(Boundaries{jj}){5} - P(end-1,:,:);
                    elseif BC.(Boundaries{jj}){5} < 0
                        % Impermeable boundary
                        P(end,:,:)    = P(end-1,:,:);
                    else
                        % Don't do anything
                    end
                    % Check isothermal boundary
                    if BC.(Boundaries{jj}){6} > 0
                        T(end,:,:)        = 2*BC.(Boundaries{jj}){6} - T(end-1,:,:);
                    else
                        % No effect of temperature
                    end
                    % Update E and rho based on P and T
                    rho  = Calculate_Rho_from_TandP( bSolver, T,P, Fluid,Substance );
                    e    = Calculate_e_from_TandPandRho( bSolver, T,rho,P, Fluid, Substance);
                    ke   = 0.5*(u.^2 + v.^2 + w.^2);
                    E    = e + ke;
                    % Update conserved variables
                    rhou = rho.*u;
                    rhov = rho.*v;
                    rhow = rho.*w;
                    rhoE = rho.*E;

                otherwise

            end

        case 'west'
            % Check type of boundary
            switch BC.(Boundaries{jj}){1}

                case 'Periodic'
                    % South boundary points
                    rho(:,1,:)     = rho(:,end-1,:);
                    rhou(:,1,:)    = rhou(:,end-1,:);
                    rhov(:,1,:)    = rhov(:,end-1,:);
                    rhow(:,1,:)    = rhow(:,end-1,:);
                    rhoE(:,1,:)    = rhoE(:,end-1,:);
                    P(:,1,:)       = P(:,end-1,:);

                case 'Neumann'
                    % South boundary points
                    rho(:,1,:)     = rho(:,2,:);
                    rhou(:,1,:)    = rhou(:,2,:);
                    rhov(:,1,:)    = rhov(:,2,:);
                    rhow(:,1,:)    = rhow(:,2,:);
                    rhoE(:,1,:)    = rhoE(:,2,:);
                    P(:,1,:)       = P(:,2,:);

                case 'Dirichlet'
                    % South boundary points
                    u(:,1,:)       = 2*BC.(Boundaries{jj}){2} - u(:,2,:);
                    v(:,1,:)       = 2*BC.(Boundaries{jj}){3} - v(:,2,:);
                    w(:,1,:)       = 2*BC.(Boundaries{jj}){4} - w(:,2,:);
                    % Check impermeable boundary
                    if BC.(Boundaries{jj}){5} > 0
                        P(:,1,:)   = 2*BC.(Boundaries{jj}){5} - P(:,2,:);
                    elseif BC.(Boundaries{jj}){5} < 0
                        % Impermeable boundary
                        P(:,1,:)    = P(:,2,:);
                    else
                        % Don't do anything
                    end
                    % Check isothermal boundary
                    if BC.(Boundaries{jj}){6} > 0
                        T(:,1,:)        = 2*BC.(Boundaries{jj}){6} - T(:,2,:);
                    else
                        % No effect of temperature
                    end
                    % Update E and rho based on P and T
                    rho  = Calculate_Rho_from_TandP( bSolver, T,P, Fluid,Substance );
                    e    = Calculate_e_from_TandPandRho( bSolver, T,rho,P, Fluid, Substance);
                    ke   = 0.5*(u.^2 + v.^2 + w.^2);
                    E    = e + ke;
                    % Update conserved variables
                    rhou = rho.*u;
                    rhov = rho.*v;
                    rhow = rho.*w;
                    rhoE = rho.*E;

                otherwise

            end

        case 'east'
            % Check type of boundary
            switch BC.(Boundaries{jj}){1}

                case 'Periodic'
                    % North boundary points
                    rho(:,end,:)   = rho(:,2,:);
                    rhou(:,end,:)  = rhou(:,2,:);
                    rhov(:,end,:)  = rhov(:,2,:);
                    rhow(:,end,:)  = rhow(:,2,:);
                    rhoE(:,end,:)  = rhoE(:,2,:);
                    P(:,end,:)     = P(:,2,:);

                case 'Neumann'
                    % North boundary points
                    rho(:,end,:)   = rho(:,end-1,:);
                    rhou(:,end,:)  = rhou(:,end-1,:);
                    rhov(:,end,:)  = rhov(:,end-1,:);
                    rhow(:,end,:)  = rhow(:,end-1,:);
                    rhoE(:,end,:)  = rhoE(:,end-1,:);
                    P(:,end,:)     = P(:,end-1,:);

                case 'Dirichlet'
                    % North boundary points
                    u(:,end,:)       = 2*BC.(Boundaries{jj}){2} - u(:,end-1,:);
                    v(:,end,:)       = 2*BC.(Boundaries{jj}){3} - v(:,end-1,:);
                    w(:,end,:)       = 2*BC.(Boundaries{jj}){4} - w(:,end-1,:);
                    % Check impermeable boundary
                    if BC.(Boundaries{jj}){5} > 0
                        P(:,end,:)   = 2*BC.(Boundaries{jj}){5} - P(:,end-1,:);
                    elseif BC.(Boundaries{jj}){5} < 0
                        % Impermeable boundary
                        P(:,end,:)    = P(:,end-1,:);
                    else
                        % Don't do anything
                    end
                    % Check isothermal boundary
                    if BC.(Boundaries{jj}){6} > 0
                        T(:,end,:)        = 2*BC.(Boundaries{jj}){6} - T(:,end-1,:);
                    else
                        % No effect of temperature
                    end
                    % Update E and rho based on P and T
                    rho  = Calculate_Rho_from_TandP( bSolver, T,P, Fluid,Substance );
                    e    = Calculate_e_from_TandPandRho( bSolver, T,rho,P, Fluid, Substance);
                    ke   = 0.5*(u.^2 + v.^2 + w.^2);
                    E    = e + ke;
                    % Update conserved variables
                    rhou = rho.*u;
                    rhov = rho.*v;
                    rhow = rho.*w;
                    rhoE = rho.*E;

                otherwise

            end

        case 'front'
            % Check type of boundary
            switch BC.(Boundaries{jj}){1}

                case 'Periodic'
                    % Back boundary points
                    rho(:,:,1)     = rho(:,:,end-1);
                    rhou(:,:,1)    = rhou(:,:,end-1);
                    rhov(:,:,1)    = rhov(:,:,end-1);
                    rhow(:,:,1)    = rhow(:,:,end-1);
                    rhoE(:,:,1)    = rhoE(:,:,end-1);
                    P(:,:,1)       = P(:,:,end-1);

                case 'Neumann'
                    % Back boundary points
                    rho(:,:,1)     = rho(:,:,2);
                    rhou(:,:,1)    = rhou(:,:,2);
                    rhov(:,:,1)    = rhov(:,:,2);
                    rhow(:,:,1)    = rhow(:,:,2);
                    rhoE(:,:,1)    = rhoE(:,:,2);
                    P(:,:,1)       = P(:,:,2);

                case 'Dirichlet'
                    % Back boundary points
                    u(:,:,1)       = 2*BC.(Boundaries{jj}){2} - u(:,:,2);
                    v(:,:,1)       = 2*BC.(Boundaries{jj}){3} - v(:,:,2);
                    w(:,:,1)       = 2*BC.(Boundaries{jj}){4} - w(:,:,2);
                    % Check impermeable boundary
                    if BC.(Boundaries{jj}){5} > 0
                        P(:,:,1)   = 2*BC.(Boundaries{jj}){5} - P(:,:,2);
                    elseif BC.(Boundaries{jj}){5} < 0
                        % Impermeable boundary
                        P(:,:,1)    = P(:,2);
                    else
                        % Don't do anything
                    end
                    % Check isothermal boundary
                    if BC.(Boundaries{jj}){6} > 0
                        T(:,:,1)        = 2*BC.(Boundaries{jj}){6} - T(:,:,2);
                    else
                        % No effect of temperature
                    end
                    % Update E and rho based on P and T
                    rho  = Calculate_Rho_from_TandP( bSolver, T,P, Fluid,Substance );
                    e    = Calculate_e_from_TandPandRho( bSolver, T,rho,P, Fluid, Substance);
                    ke   = 0.5*(u.^2 + v.^2 + w.^2);
                    E    = e + ke;
                    % Update conserved variables
                    rhou = rho.*u;
                    rhov = rho.*v;
                    rhow = rho.*w;
                    rhoE = rho.*E;

                otherwise

            end

        case 'back'
            % Check type of boundary
            switch BC.(Boundaries{jj}){1}

                case 'Periodic'
                    % Front boundary points
                    rho(:,:,end)     = rho(:,:,2);
                    rhou(:,:,end)    = rhou(:,:,2);
                    rhov(:,:,end)    = rhov(:,:,2);
                    rhow(:,:,end)    = rhow(:,:,2);
                    rhoE(:,:,end)    = rhoE(:,:,2);
                    P(:,:,end)       = P(:,:,2);

                case 'Neumann'
                    % Front boundary points
                    rho(:,:,end)   = rho(:,:,end-1);
                    rhou(:,:,end)  = rhou(:,:,end-1);
                    rhov(:,:,end)  = rhov(:,:,end-1);
                    rhow(:,:,end)  = rhow(:,:,end-1);
                    rhoE(:,:,end)  = rhoE(:,:,end-1);
                    P(:,:,end)     = P(:,:,end-1);

                case 'Dirichlet'
                    % South boundary points
                    u(:,:,end)       = 2*BC.(Boundaries{jj}){2} - u(:,:,end-1);
                    v(:,:,end)       = 2*BC.(Boundaries{jj}){3} - v(:,:,end-1);
                    w(:,:,end)       = 2*BC.(Boundaries{jj}){4} - w(:,:,end-1);
                    % Check impermeable boundary
                    if BC.(Boundaries{jj}){5} > 0
                        P(:,:,end)   = 2*BC.(Boundaries{jj}){5} - P(:,:,end-1);
                    elseif BC.(Boundaries{jj}){5} < 0
                        % Impermeable boundary
                        P(:,:,end)    = P(:,:,end-1);
                    else
                        % Don't do anything
                    end
                    % Check isothermal boundary
                    if BC.(Boundaries{jj}){6} > 0
                        T(:,:,end)        = 2*BC.(Boundaries{jj}){6} - T(:,:,end-1);
                    else
                        % No effect of temperature
                    end
                    % Update E and rho based on P and T
                    rho  = Calculate_Rho_from_TandP( bSolver, T,P, Fluid,Substance );
                    e    = Calculate_e_from_TandPandRho( bSolver, T,rho,P, Fluid, Substance);
                    ke   = 0.5*(u.^2 + v.^2 + w.^2);
                    E    = e + ke;

                    % Update conserved variables
                    rhou = rho.*u;
                    rhov = rho.*v;
                    rhow = rho.*w;
                    rhoE = rho.*E;

                otherwise

            end

        otherwise


    end

end

%% Update axis
% Y axis
j = 1; k = 1;
rho(2:end-1,j,k)  = 0.5*(rho(2:end-1,j+1,k)  + rho(2:end-1,j,k+1));
rhou(2:end-1,j,k) = 0.5*(rhou(2:end-1,j+1,k) + rhou(2:end-1,j,k+1));
rhov(2:end-1,j,k) = 0.5*(rhov(2:end-1,j+1,k) + rhov(2:end-1,j,k+1));
rhow(2:end-1,j,k) = 0.5*(rhow(2:end-1,j+1,k) + rhow(2:end-1,j,k+1));
rhoE(2:end-1,j,k) = 0.5*(rhoE(2:end-1,j+1,k) + rhoE(2:end-1,j,k+1));
P(2:end-1,j,k)    = 0.5*(P(2:end-1,j+1,k) + P(2:end-1,j,k+1));

j = 1; k = length(rho(1,1,:));
rho(2:end-1,j,k)  = 0.5*(rho(2:end-1,j+1,k)  + rho(2:end-1,j,k-1));
rhou(2:end-1,j,k) = 0.5*(rhou(2:end-1,j+1,k) + rhou(2:end-1,j,k-1));
rhov(2:end-1,j,k) = 0.5*(rhov(2:end-1,j+1,k) + rhov(2:end-1,j,k-1));
rhow(2:end-1,j,k) = 0.5*(rhow(2:end-1,j+1,k) + rhow(2:end-1,j,k-1));
rhoE(2:end-1,j,k) = 0.5*(rhoE(2:end-1,j+1,k) + rhoE(2:end-1,j,k-1));
P(2:end-1,j,k)    = 0.5*(P(2:end-1,j+1,k) + P(2:end-1,j,k-1));

j = length(rho(1,:,1)); k = 1;
rho(2:end-1,j,k)  = 0.5*(rho(2:end-1,j-1,k)  + rho(2:end-1,j,k+1));
rhou(2:end-1,j,k) = 0.5*(rhou(2:end-1,j-1,k) + rhou(2:end-1,j,k+1));
rhov(2:end-1,j,k) = 0.5*(rhov(2:end-1,j-1,k) + rhov(2:end-1,j,k+1));
rhow(2:end-1,j,k) = 0.5*(rhow(2:end-1,j-1,k) + rhow(2:end-1,j,k+1));
rhoE(2:end-1,j,k) = 0.5*(rhoE(2:end-1,j-1,k) + rhoE(2:end-1,j,k+1));
P(2:end-1,j,k)    = 0.5*(P(2:end-1,j-1,k) + P(2:end-1,j,k+1));

j = length(rho(1,:,1)); k = length(rho(1,1,:));
rho(2:end-1,j,k)  = 0.5*(rho(2:end-1,j-1,k)  + rho(2:end-1,j,k-1));
rhou(2:end-1,j,k) = 0.5*(rhou(2:end-1,j-1,k) + rhou(2:end-1,j,k-1));
rhov(2:end-1,j,k) = 0.5*(rhov(2:end-1,j-1,k) + rhov(2:end-1,j,k-1));
rhow(2:end-1,j,k) = 0.5*(rhow(2:end-1,j-1,k) + rhow(2:end-1,j,k-1));
rhoE(2:end-1,j,k) = 0.5*(rhoE(2:end-1,j-1,k) + rhoE(2:end-1,j,k-1));
P(2:end-1,j,k)    = 0.5*(P(2:end-1,j-1,k) + P(2:end-1,j,k-1));

% X axis
i = 1; k = 1;
rho(i,2:end-1,k)  = 0.5*(rho(i+1,2:end-1,k)  + rho(i,2:end-1,k+1));
rhou(i,2:end-1,k) = 0.5*(rhou(i+1,2:end-1,k) + rhou(i,2:end-1,k+1));
rhov(i,2:end-1,k) = 0.5*(rhov(i+1,2:end-1,k) + rhov(i,2:end-1,k+1));
rhow(i,2:end-1,k) = 0.5*(rhow(i+1,2:end-1,k) + rhow(i,2:end-1,k+1));
rhoE(i,2:end-1,k) = 0.5*(rhoE(i+1,2:end-1,k) + rhoE(i,2:end-1,k+1));
P(i,2:end-1,k)    = 0.5*(P(i+1,2:end-1,k) + P(i,2:end-1,k+1));

i = 1; k = length(rho(1,1,:));
rho(i,2:end-1,k)  = 0.5*(rho(i+1,2:end-1,k)  + rho(i,2:end-1,k-1));
rhou(i,2:end-1,k) = 0.5*(rhou(i+1,2:end-1,k) + rhou(i,2:end-1,k-1));
rhov(i,2:end-1,k) = 0.5*(rhov(i+1,2:end-1,k) + rhov(i,2:end-1,k-1));
rhow(i,2:end-1,k) = 0.5*(rhow(i+1,2:end-1,k) + rhow(i,2:end-1,k-1));
rhoE(i,2:end-1,k) = 0.5*(rhoE(i+1,2:end-1,k) + rhoE(i,2:end-1,k-1));
P(i,2:end-1,k)    = 0.5*(P(i+1,2:end-1,k) + P(i,2:end-1,k-1));

i = length(rho(:,1,1)); k = 1;
rho(i,2:end-1,k)  = 0.5*(rho(i-1,2:end-1,k)  + rho(i,2:end-1,k+1));
rhou(i,2:end-1,k) = 0.5*(rhou(i-1,2:end-1,k) + rhou(i,2:end-1,k+1));
rhov(i,2:end-1,k) = 0.5*(rhov(i-1,2:end-1,k) + rhov(i,2:end-1,k+1));
rhow(i,2:end-1,k) = 0.5*(rhow(i-1,2:end-1,k) + rhow(i,2:end-1,k+1));
rhoE(i,2:end-1,k) = 0.5*(rhoE(i-1,2:end-1,k) + rhoE(i,2:end-1,k+1));
P(i,2:end-1,k)    = 0.5*(P(i-1,2:end-1,k) + P(i,2:end-1,k+1));

i = length(rho(:,1,1)); k = length(rho(1,1,:));
rho(i,2:end-1,k)  = 0.5*(rho(i-1,2:end-1,k)  + rho(i,2:end-1,k-1));
rhou(i,2:end-1,k) = 0.5*(rhou(i-1,2:end-1,k) + rhou(i,2:end-1,k-1));
rhov(i,2:end-1,k) = 0.5*(rhov(i-1,2:end-1,k) + rhov(i,2:end-1,k-1));
rhow(i,2:end-1,k) = 0.5*(rhow(i-1,2:end-1,k) + rhow(i,2:end-1,k-1));
rhoE(i,2:end-1,k) = 0.5*(rhoE(i-1,2:end-1,k) + rhoE(i,2:end-1,k-1));
P(i,2:end-1,k)    = 0.5*(P(i-1,2:end-1,k) + P(i,2:end-1,k-1));

% Z axis
i = 1; j = 1;
rho(i,j,2:end-1)  = 0.5*(rho(i+1,j,2:end-1)  + rho(i,j+1,2:end-1));
rhou(i,j,2:end-1) = 0.5*(rhou(i+1,j,2:end-1) + rhou(i,j+1,2:end-1));
rhov(i,j,2:end-1) = 0.5*(rhov(i+1,j,2:end-1) + rhov(i,j+1,2:end-1));
rhow(i,j,2:end-1) = 0.5*(rhow(i+1,j,2:end-1) + rhow(i,j+1,2:end-1));
rhoE(i,j,2:end-1) = 0.5*(rhoE(i+1,j,2:end-1) + rhoE(i,j+1,2:end-1));
P(i,j,2:end-1)    = 0.5*(P(i+1,j,2:end-1) + P(i,j+1,2:end-1));

i = 1; j = length(rho(1,:,1));
rho(i,j,2:end-1)  = 0.5*(rho(i+1,j,2:end-1)  + rho(i,j-1,2:end-1));
rhou(i,j,2:end-1) = 0.5*(rhou(i+1,j,2:end-1) + rhou(i,j-1,2:end-1));
rhov(i,j,2:end-1) = 0.5*(rhov(i+1,j,2:end-1) + rhov(i,j-1,2:end-1));
rhow(i,j,2:end-1) = 0.5*(rhow(i+1,j,2:end-1) + rhow(i,j-1,2:end-1));
rhoE(i,j,2:end-1) = 0.5*(rhoE(i+1,j,2:end-1) + rhoE(i,j-1,2:end-1));
P(i,j,2:end-1)    = 0.5*(P(i+1,j,2:end-1) + P(i,j-1,2:end-1));

i = length(rho(:,1,1)); j = 1;
rho(i,j,2:end-1)  = 0.5*(rho(i-1,j,2:end-1)  + rho(i,j+1,2:end-1));
rhou(i,j,2:end-1) = 0.5*(rhou(i-1,j,2:end-1) + rhou(i,j+1,2:end-1));
rhov(i,j,2:end-1) = 0.5*(rhov(i-1,j,2:end-1) + rhov(i,j+1,2:end-1));
rhow(i,j,2:end-1) = 0.5*(rhow(i-1,j,2:end-1) + rhow(i,j+1,2:end-1));
rhoE(i,j,2:end-1) = 0.5*(rhoE(i-1,j,2:end-1) + rhoE(i,j+1,2:end-1));
P(i,j,2:end-1)    = 0.5*(P(i-1,j,2:end-1) + P(i,j+1,2:end-1));

i = length(rho(:,1,1)); j = length(rho(1,:,1));
rho(i,j,2:end-1)  = 0.5*(rho(i-1,j,2:end-1)  + rho(i,j-1,2:end-1));
rhou(i,j,2:end-1) = 0.5*(rhou(i-1,j,2:end-1) + rhou(i,j-1,2:end-1));
rhov(i,j,2:end-1) = 0.5*(rhov(i-1,j,2:end-1) + rhov(i,j-1,2:end-1));
rhow(i,j,2:end-1) = 0.5*(rhow(i-1,j,2:end-1) + rhow(i,j-1,2:end-1));
rhoE(i,j,2:end-1) = 0.5*(rhoE(i-1,j,2:end-1) + rhoE(i,j-1,2:end-1));
P(i,j,2:end-1)    = 0.5*(P(i-1,j,2:end-1) + P(i,j-1,2:end-1));

%% Corners
% i = 1, j = 1, k = 1
rho(1,1,1)  = 1/3*(rho(2,1,1)  + rho(1,2,1)  + rho(1,1,2));
rhou(1,1,1) = 1/3*(rhou(2,1,1) + rhou(1,2,1) + rhou(1,1,2));
rhov(1,1,1) = 1/3*(rhov(2,1,1) + rhov(1,2,1) + rhov(1,1,2));
rhow(1,1,1) = 1/3*(rhow(2,1,1) + rhow(1,2,1) + rhow(1,1,2));
rhoE(1,1,1) = 1/3*(rhoE(2,1,1) + rhoE(1,2,1) + rhoE(1,1,2));
P(1,1,1)    = 1/3*(P(2,1,1) + P(1,2,1) + P(1,1,2));

% i = end, j = 1, k = 1
rho(end,1,1)  = 1/3*(rho(end-1,1,1)  + rho(end,2,1)  + rho(end,1,2));
rhou(end,1,1) = 1/3*(rhou(end-1,1,1) + rhou(end,2,1) + rhou(end,1,2));
rhov(end,1,1) = 1/3*(rhov(end-1,1,1) + rhov(end,2,1) + rhov(end,1,2));
rhow(end,1,1) = 1/3*(rhow(end-1,1,1) + rhow(end,2,1) + rhow(end,1,2));
rhoE(end,1,1) = 1/3*(rhoE(end-1,1,1) + rhoE(end,2,1) + rhoE(end,1,2));
P(end,1,1)    = 1/3*(P(end-1,1,1) + P(end,2,1) + P(end,1,2));

% i = 1, j = end, k = 1
rho(1,end,1)  = 1/3*(rho(2,end,1)  + rho(1,end-1,1)  + rho(1,end,2));
rhou(1,end,1) = 1/3*(rhou(2,end,1) + rhou(1,end-1,1) + rhou(1,end,2));
rhov(1,end,1) = 1/3*(rhov(2,end,1) + rhov(1,end-1,1) + rhov(1,end,2));
rhow(1,end,1) = 1/3*(rhow(2,end,1) + rhow(1,end-1,1) + rhow(1,end,2));
rhoE(1,end,1) = 1/3*(rhoE(2,end,1) + rhoE(1,end-1,1) + rhoE(1,end,2));
P(1,end,1)    = 1/3*(P(2,end,1) + P(1,end-1,1) + P(1,end,2));


% i = end, j = end, k = 1
rho(end,end,1)  = 1/3*(rho(end-1,end,1)  + rho(end,end-1,1)  + rho(end,end,2));
rhou(end,end,1) = 1/3*(rhou(end-1,end,1) + rhou(end,end-1,1) + rhou(end,end,2));
rhov(end,end,1) = 1/3*(rhov(end-1,end,1) + rhov(end,end-1,1) + rhov(end,end,2));
rhow(end,end,1) = 1/3*(rhow(end-1,end,1) + rhow(end,end-1,1) + rhow(end,end,2));
rhoE(end,end,1) = 1/3*(rhoE(end-1,end,1) + rhoE(end,end-1,1) + rhoE(end,end,2));
P(end,end,1)    = 1/3*(P(end-1,end,1) + P(end,end-1,1) + P(end,end,2));

% i = 1, j = 1, k = end
rho(1,1,end)  = 1/3*(rho(2,1,end)  + rho(1,2,end)  + rho(1,1,end-1));
rhou(1,1,end) = 1/3*(rhou(2,1,end) + rhou(1,2,end) + rhou(1,1,end-1));
rhov(1,1,end) = 1/3*(rhov(2,1,end) + rhov(1,2,end) + rhov(1,1,end-1));
rhow(1,1,end) = 1/3*(rhow(2,1,end) + rhow(1,2,end) + rhow(1,1,end-1));
rhoE(1,1,end) = 1/3*(rhoE(2,1,end) + rhoE(1,2,end) + rhoE(1,1,end-1));
P(1,1,end)    = 1/3*(P(2,1,end) + P(1,2,end) + P(1,1,end-1));

% i = end, j = 1, k = end
rho(end,1,end)  = 1/3*(rho(end-1,1,end)  + rho(end,2,end)  + rho(end,1,end-1));
rhou(end,1,end) = 1/3*(rhou(end-1,1,end) + rhou(end,2,end) + rhou(end,1,end-1));
rhov(end,1,end) = 1/3*(rhov(end-1,1,end) + rhov(end,2,end) + rhov(end,1,end-1));
rhow(end,1,end) = 1/3*(rhow(end-1,1,end) + rhow(end,2,end) + rhow(end,1,end-1));
rhoE(end,1,end) = 1/3*(rhoE(end-1,1,end) + rhoE(end,2,end) + rhoE(end,1,end-1));
P(end,1,end)    = 1/3*(P(end-1,1,end) + P(end,2,end) + P(end,1,end-1));

% i = 1, j = end, k = end
rho(1,end,end)  = 1/3*(rho(2,end,end)  + rho(1,end-1,end)  + rho(1,end,end-1));
rhou(1,end,end) = 1/3*(rhou(2,end,end) + rhou(1,end-1,end) + rhou(1,end,end-1));
rhov(1,end,end) = 1/3*(rhov(2,end,end) + rhov(1,end-1,end) + rhov(1,end,end-1));
rhow(1,end,end) = 1/3*(rhow(2,end,end) + rhow(1,end-1,end) + rhow(1,end,end-1));
rhoE(1,end,end) = 1/3*(rhoE(2,end,end) + rhoE(1,end-1,end) + rhoE(1,end,end-1));
P(1,end,end)    = 1/3*(P(2,end,end) + P(1,end-1,end) + P(1,end,end-1));

% i = end, j = end, k = end
rho(end,end,end)  = 1/3*(rho(end-1,end,end)  + rho(end,end-1,end)  + rho(end,end,end-1));
rhou(end,end,end) = 1/3*(rhou(end-1,end,end) + rhou(end,end-1,end) + rhou(end,end,end-1));
rhov(end,end,end) = 1/3*(rhov(end-1,end,end) + rhov(end,end-1,end) + rhov(end,end,end-1));
rhow(end,end,end) = 1/3*(rhow(end-1,end,end) + rhow(end,end-1,end) + rhow(end,end,end-1));
rhoE(end,end,end) = 1/3*(rhoE(end-1,end,end) + rhoE(end,end-1,end) + rhoE(end,end,end-1));
P(end,end,end)    = 1/3*(P(end-1,end,end) + P(end,end-1,end) + P(end,end,end-1));


end

