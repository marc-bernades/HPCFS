function [u,v,w,E,rho,mu,kappa,c_v,c_p,P,T,ke,e,sos,Beta_T, Beta_v, Beta_s, Alpha_p,time,t_vec,ke_total,X,Y,Z, Invariants] = Compressible_Solver(Fluid,bSolver,x_0,y_0,z_0,L_x,L_y,L_z,t_0,t_end,name_file_out,...
    num_grid_x,num_grid_y,num_grid_z, CFL, n_iter_max, output_iter, RK_order, scheme, BC, ST, K, bViscous_5PointsStencil, bViscosityVariable, bKappaVariable, bEnergySplit, bEnergySplitEnthalpy, Test, Substance, HP_model, bPressureModel, bFlux, bFilter, varargin)

%format long

%% VARIABLE INITIALIZATION
[grid_mesh, rho, u, v, w, E, rhou, rhov, rhow, rhoE, P, T, mu, kappa, sos, rhou_0, rhov_0, rhow_0, rhoE_0, ...
    rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv, dP_rhou, dP_rhov, dP_rhow, dP_rhoE, ...
    rhou_visc, rhov_visc, rhow_visc, rhoE_visc, f_rhou, f_rhov, f_rhow, f_rhoE] = Initialize_Variables(num_grid_x, num_grid_y, num_grid_z);

%% CREATE GRID
% Uniform grid
grid_x = (x_0 + -0.5*L_x/num_grid_x):L_x/num_grid_x:(L_x + x_0 + 0.5*L_x/num_grid_x);         % Vector of x grid points [m]
grid_y = (y_0 + -0.5*L_y/num_grid_y):L_y/num_grid_y:(L_y + y_0 + 0.5*L_y/num_grid_y);         % Vector of y grid points [m]
grid_z = (z_0 + -0.5*L_z/num_grid_z):L_z/num_grid_z:(L_z + z_0 + 0.5*L_z/num_grid_z);         % Vector of z grid points [m]

% Stretching grid
try
    if Fluid.A_x ~= 0
        grid_x = stretching_grid(num_grid_x,x_0,Fluid.A_x,L_x);
        disp('Stretching discretization in X-direction...')
    end
    if Fluid.A_y ~= 0
        grid_y = stretching_grid(num_grid_y,y_0,Fluid.A_y,L_y);
        disp('Stretching discretization in Y-direction...')
    end
    if Fluid.A_z ~= 0
        grid_z = stretching_grid(num_grid_z,z_0,Fluid.A_z,L_z);
        disp('Stretching discretization in Z-direction...')
    end
catch
    disp('Uniform grid discretization...')
end

% Mesh
[X,Y,Z] = meshgrid(grid_x,grid_y,grid_z);
N = max([num_grid_x,num_grid_y,num_grid_z]);
L = max([L_x,L_y,L_z]);

% Delta x, y and z based on CentralDerivative_d1_2ndOrder
[dx,~,~] = CentralDerivative_d1_2ndOrder(X);
[~,dy,~] = CentralDerivative_d1_2ndOrder(Y);
[~,~,dz] = CentralDerivative_d1_2ndOrder(Z);

%% Initialization
% Initialize velocity field and P and T
[u,v,w,P,T] = Initialize_Fields(bSolver,X,Y,Z,L_x,L_y,L_z,Test,Fluid,Substance);

%Initialize thermodynamics
[rho,~,~,E,~,~,~,~,~,~,~,~]  = Initialize_Thermodynamics(bSolver, P,T,u,v,w,bPressureModel,Fluid,Substance);

% Update thermodynamic state
[ke, e, ~, ~, c_p, c_v, gamma, sos, Beta_T, Beta_v, Beta_s, Alpha_p,~] = Update_ThermodynamicState(bSolver,rho,u,v,w,E,P,T,Fluid,Substance,bPressureModel);

% Initialize viscous terms
[mu,kappa]                   = Calculate_HighPressure_Transport_Coeff(bSolver, mu, kappa, Fluid, Substance, T, rho, HP_model);

% Update source terms based on initial conditions
f_rhou   = f_rhou + ST.f_rhou;
f_rhov   = f_rhov + ST.f_rhov;
f_rhow   = f_rhow + ST.f_rhow;
f_rhoE   = f_rhoE + ST.f_rhoE;

%Initialize time
time = t_0;

% Filter matrix
%[A,B] = Filter_Matrix_F2(u);
[A,B] = Filter_Matrix_F4(u);


% Check Boundary Conditions type to feed sweep for Dirichlet T and P
% preparation
Boundaries = fieldnames(BC);
BC_Type    = cell(1,length(Boundaries));
for jj = 1:length(Boundaries)
    BC_Type{jj} = BC.(Boundaries{jj}){1};
end

%% Update conserved and terms
% Update conserved variables
% [rhou, rhov, rhow, rhoE] = Update_ConservedVariables(rho,u,v,w,E);
[rhou, rhov, rhow, rhoE] = Update_ConservedVariables(rho,u,v,w,e);


% Update current state
[rho_0, rhou_0, rhov_0, rhow_0, rhoE_0] = Update_CurrentState(rho, rhou, rhov, rhow, rhoE);
P_0  = P;
ke_0 = sum(sum(sum(ke(2:end-1,2:end-1,2:end-1)))); % For print

% Pre-allocate memory Initialize time variation variables and time array
if isempty(varargin)
    dt = Calculate_TimeStep(bSolver,rho,rhou,rhov,rhow,CFL,sos,gamma,c_p,kappa,mu,N, L,dx,dy,dz,Fluid,Test);
else
    dt = varargin{1};
end

% Initialize Invariants and time dependence quantities
b_flag   = 11;
margin   = 100;
[t_vec, ke_total, Invariants] = Initialize_Invariants(b_flag,margin,t_end,t_0,dt);
Ref      = [];



%% Time integrator coefficients
[a,b,c] = ButcherTableau(RK_order);

%% Final time step to obtain exact t_end
dt_end = 0;
tic
try
    % Time iteration
    for t = 1:n_iter_max
        
        % Time evolution of linear and quadratic invariants
        t_vec(t)         = t_vec(t) + time;
        ke_total(t)      = Calculate_ke(u, v, w, X, Y, Z);
        [Invariants,Ref] = Calculate_Invariants(dx,dy,dz,rho, u, v, w, E, e, c_v , gamma, T, P,t,Invariants,Ref,Substance,Test,bSolver);

        % Calculate time step
        if dt_end == 0
            if isempty(varargin)
                dt = Calculate_TimeStep(bSolver,rho,rhou,rhov,rhow,CFL,sos,gamma,c_p,kappa,mu,N, L,dx,dy,dz,Fluid,Test);
                disp("Computing iteration " + num2str(t) + " at t = " + num2str(time) + "s" + " dt = " + num2str(dt) + " ke/ke0 = " + num2str(sum(sum(sum(ke(2:end-1,2:end-1,2:end-1))))/ke_0))
            else
                dt = varargin{1};
                U_ref = sqrt(max(max(max((rhou.*rhou + rhov.*rhov + rhow.*rhow)./rho.^2)))) + max(max(max(sos)));
                CFL = dt/((L/N)/U_ref);
                disp("Computing iteration " + num2str(t) + " at t = " + num2str(time) + "s" + " dt = " + num2str(dt) + " and CFL = " + num2str(CFL) + "ke/ke0 = " + num2str(sum(sum(sum(ke(2:end-1,2:end-1,2:end-1))))/ke_0))
            end
        else
            % Adjust final time step to match t_end
            dt = dt_end;
        end

        % Time integration - RK

        for i_RK = 1:RK_order
            % Obtain source terms

            % Calculate convective terms
            % Matricial > Energy split e + ke + p/rho
            if bEnergySplit == 1
                % General method
%                 [C_D, C_u, C_fi, C_rho, C_L] = ConvectiveTermSplit(rho,u,v,w,E,X,Y,Z,P); %Swap E to e for PEP-RG
                [C_D, C_u, C_fi, C_rho, C_L] = ConvectiveTermSplit(rho,u,v,w,e,X,Y,Z,P); %Swap E to e for PEP-RG
                % Calculate convective terms based on coefficients defined
                rho_conv  = Calculate_Convective_Split(C_D{1}, C_u{1}, C_fi{1}, C_rho{1}, C_L{1}, scheme.mass);
                rhou_conv = Calculate_Convective_Split(C_D{2}, C_u{2}, C_fi{2}, C_rho{2}, C_L{2}, scheme.momentum);
                rhov_conv = Calculate_Convective_Split(C_D{3}, C_u{3}, C_fi{3}, C_rho{3}, C_L{3}, scheme.momentum);
                rhow_conv = Calculate_Convective_Split(C_D{4}, C_u{4}, C_fi{4}, C_rho{4}, C_L{4}, scheme.momentum);
                rhoe_conv = Calculate_Convective_Split(C_D{5}, C_u{5}, C_fi{5}, C_rho{5}, C_L{5}, scheme.e);
                rhok_conv = Calculate_Convective_Split(C_D{6}, C_u{6}, C_fi{6}, C_rho{6}, C_L{6}, scheme.k);
                Pu_conv   = Calculate_Convective_Split(C_D{7}, C_u{7}, C_fi{7}, C_rho{7}, C_L{7}, scheme.p);
                rhoE_conv = rhoe_conv + rhok_conv + Pu_conv;

                % Matricial method specific
                if bEnergySplitEnthalpy == 1
                    % Enthalpy based including pressure gradient and energy
                    [C_D, C_u, C_fi, C_rho, C_L] = ConvectiveTermEnthalpySplit(rho,u,v,w,E,rhou,rhov,rhow,rhoE,X,Y,Z,P);
                    % Calculate convective terms based on scheme defined
                    [rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv] = Calculate_Convective(C_D, C_u, C_fi, C_rho, C_L, scheme.def);

                elseif bEnergySplitEnthalpy == 0 && bEnergySplit == 0
                    % Energy equation (standard NS equations)
                    [C_D, C_u, C_fi, C_rho, C_L] = ConvectiveTerm(rho,u,v,w,E,rhou,rhov,rhow,rhoE,X,Y,Z);
                    % Calculate convective terms based on scheme defined
                    [rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv] = Calculate_Convective(C_D, C_u, C_fi, C_rho, C_L, scheme.def);
                end
            end

            % Calculate pressure gradient
            [dP_rhou, dP_rhov, dP_rhow, dP_rhoE] = PressureGradient(u,v,w,P,dx,dy,dz);

            % Activate sensor if needed
            if bFlux.hybrid == 1
                % Apply filter to the energy split convective term
                fi = Sensor_Inviscid(u,v,w,dx,dy,dz);
                drawnow
                plot(fi(2,:,2)) 
                hold on
                rho_conv_hyb = (1-fi).*rho_conv;
                rhou_conv_hyb = (1-fi).*(rhou_conv  + dP_rhou);
                rhov_conv_hyb = (1-fi).*(rhov_conv  + dP_rhov);
                rhow_conv_hyb = (1-fi).*(rhow_conv  + dP_rhow);

                if bEnergySplit == 1 ||  bEnergySplitEnthalpy == 1
                    % rhoE includes pressure gradient
                    rhoE_conv_hyb = (1-fi).*rhoE_conv;
                else
                    rhoE_conv_hyb = (1-fi).*(rhoE_conv + dP_rhoE);
                end


            end
                % Flux method
            if bFlux.gate == 1
                if strcmp(bFlux.type,'WENO5')
                    % WENO5 flux form
                    [rho_conv,rhou_conv,rhov_conv,rhow_conv,rhoE_conv] = CalculateInviscidFlux_WENO5(rho,u,v,w,E,sos,P,dx,dy,dz,bFlux, bSolver, Fluid, Substance);
                else
                    % Flux form for KGP, Shima, HLLC, HLLC+, etc
                    [rho_conv,rhou_conv,rhov_conv,rhow_conv,rhoE_conv] = CalculateInviscidFlux(rho,u,v,w,E,sos,P,dx,dy,dz,bFlux, bSolver, Fluid, Substance);               
                end
            end

            if bFlux.hybrid == 1
                % Apply filter to the energy split convective term
                rho_conv  = rho_conv_hyb  + (fi).*rho_conv;
                rhou_conv = rhou_conv_hyb + (fi).*rhou_conv;
                rhov_conv = rhov_conv_hyb + (fi).*rhov_conv;
                rhow_conv = rhow_conv_hyb + (fi).*rhow_conv;
                rhoE_conv = rhoE_conv_hyb + (fi).*rhoE_conv;
            end

            %For 1D_Adv we need to decouple convective term using correctly u.
            if strcmp(Test,'1D_Adv')
                if strcmp(bSolver,'Ideal')
                    % Only needed for Ideal gas
                    [drho_x,~,~]   = CentralDerivative_d1_2ndOrder(rho);
                    [dx,~,~]       = CentralDerivative_d1_2ndOrder(X);
                    rho_conv       = drho_x./dx;
                end
            end

                       

            % Calculate viscous term
            if bViscous_5PointsStencil == 0
                % 3-Points Stencil approach (double derivative)
                [rhou_visc, rhov_visc, rhow_visc, rhoE_visc, Tau, q_div] = Calculate_Viscous(u,v,w,mu,kappa,T,dx,dy,dz,bViscosityVariable,bKappaVariable);
            else
                % 5-Points Stencil approach (gradient of gradient)
                [rhou_visc, rhov_visc, rhow_visc, rhoE_visc] = Calculate_Viscous_5PointStencil(u,v,w,mu,kappa,T,dx,dy,dz);
            end

            % Calculate source terms
            [f_rhou, f_rhov, f_rhow, f_rhoE, f_rhouvw] = CalculateSourceTerms(f_rhou, f_rhov, f_rhow, f_rhoE,u, v, w, X, Y, Z, K, Fluid);
            

            % Sum all fluxes
            [rho_tot{i_RK}, rhou_tot{i_RK}, rhov_tot{i_RK}, rhow_tot{i_RK}, rhoE_tot{i_RK}] = Calculate_TotalFlux(rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoe_conv, ...
                dP_rhou, dP_rhov, dP_rhow, dP_rhoE,...
                rhou_visc, rhov_visc, rhow_visc, rhoE_visc,  bEnergySplit, bEnergySplitEnthalpy, bFlux.gate, ...
                f_rhou, f_rhov, f_rhow, f_rhoE, f_rhouvw); %Swap rhoE_conv for rhoe_conv for PEP-RG
            % Advance conserved variable in time
            if bPressureModel == 1
                % Update pressure gradient according to Kawai
                dP{i_RK} = Pressure_variation(P,rho,u,v,w,sos,dx,dy,dz,c_v,Alpha_p,Beta_T,Tau,q_div);
                P  = TimeIntegration_RKGeneral(P_0,dt,dP,i_RK,RK_order,a,b,c);
            end
            rho  = TimeIntegration_RKGeneral(rho_0, dt,rho_tot, i_RK,RK_order,a,b,c);
            rhou = TimeIntegration_RKGeneral(rhou_0,dt,rhou_tot,i_RK,RK_order,a,b,c);
            rhov = TimeIntegration_RKGeneral(rhov_0,dt,rhov_tot,i_RK,RK_order,a,b,c);
            rhow = TimeIntegration_RKGeneral(rhow_0,dt,rhow_tot,i_RK,RK_order,a,b,c);
            rhoE = TimeIntegration_RKGeneral(rhoE_0,dt,rhoE_tot,i_RK,RK_order,a,b,c);

            
            if bFilter.Gate == 1
                [rho]  = Filter_Conserved(rho,A,B,dx,dy,dz,bFilter.Type);
                [rhou] = Filter_Conserved(rhou,A,B,dx,dy,dz,bFilter.Type);
                [rhov] = Filter_Conserved(rhov,A,B,dx,dy,dz,bFilter.Type);
                [rhow] = Filter_Conserved(rhow,A,B,dx,dy,dz,bFilter.Type);
                [rhoE] = Filter_Conserved(rhoE,A,B,dx,dy,dz,bFilter.Type);
                if bPressureModel == 1
                    [P] = Filter_Conserved(P,A,B,dx,dy,dz,bFilter.Type);
                end
            end

            % This makes velocities being 0, but as rho is decoupled so it
            % makes this redundant
            if strcmp(Test,'1D_Adv')
                if strcmp(bSolver,'Ideal')
                    rhou = rhou*0;
                    rhov = rhou*0;
                    rhow = rhou*0;
                    rhoE = rhou*0;
                    u    = u*0 + Fluid.U_0;
                end
            end

            % Check whether Dirichet BC are used to update P and T
            if max(strcmp(BC_Type, 'Dirichlet')) > 0
                % Update primative variables for Dirichlet if needed
                [u, v, w, E]    = Update_PrimativeVariables(rho,rhou,rhov,rhow,rhoE);

                % Update thermodynamics
                [ke, e, P, T, c_p, c_v, gamma, sos, Beta_T, Beta_v, Beta_s, Alpha_p,E] = Update_ThermodynamicState(bSolver,rho,u,v,w,E,P,T,Fluid,Substance,bPressureModel);

            end

      
            % Update boundaries
            [rho, rhou, rhov, rhow, rhoE,P] = Update_Boundaries(bSolver,rho, rhou, rhov, rhow, rhoE, u, v, w, P, T, BC, Fluid, Substance);

            % Update primative variables
            [u, v, w, E] = Update_PrimativeVariables(rho,rhou,rhov,rhow,rhoE);

            if strcmp(Test,'1D_Adv')
                if strcmp(bSolver,'Ideal')
                    u  = u*0 + Fluid.U_0;
                end
            end

            % Update thermodynamic state
            [ke, e, P, T, c_p, c_v, gamma, sos, Beta_T, Beta_v, Beta_s, Alpha_p,E] = Update_ThermodynamicState(bSolver,rho,u,v,w,E,P,T,Fluid,Substance,bPressureModel);

            % Update high pressure variables
            [mu,kappa] = Calculate_HighPressure_Transport_Coeff(bSolver, mu, kappa, Fluid, Substance, T, rho, HP_model);

        end

%         % Apply filter only to conserved variables before BC and Thermo
%         if bFilter.Gate == 1
%             [rho]  = Filter_Conserved(rho,A,B,dx,dy,dz,bFilter.Type);
%             [rhou] = Filter_Conserved(rhou,A,B,dx,dy,dz,bFilter.Type);
%             [rhov] = Filter_Conserved(rhov,A,B,dx,dy,dz,bFilter.Type);
%             [rhow] = Filter_Conserved(rhow,A,B,dx,dy,dz,bFilter.Type);
%             [rhoE] = Filter_Conserved(rhoE,A,B,dx,dy,dz,bFilter.Type);
%             % Update boundaries
%             [rho, rhou, rhov, rhow, rhoE,~] = Update_Boundaries(bSolver,rho, rhou, rhov, rhow, rhoE, u, v, w, P, T, BC, Fluid, Substance);
%         end


        % Update current state
        [rho_0, rhou_0, rhov_0, rhow_0, rhoE_0] = Update_CurrentState(rho, rhou, rhov, rhow, rhoE);
        if bPressureModel == 1
            P_0 = P;
        end

        % Update time
        time = time + dt;

        % Final time step to match t_end
        if abs(time-t_end)<dt
            dt_end = abs(time-t_end);
        end

        % Simulation completed
        if time >= t_end
            if ~exist('tt','var')
                tt = 1;
                t_out = output_iter(tt);
            end
            break
        end

        % Compute spectra > Output files 
        t_out = output_iter;
        for tt = 1:length(t_out)
            if abs(time - t_out(tt))<dt
                DataOutput(u,v,w,E,rho,mu,kappa,P,T,ke,e,c_p,c_v,sos,Beta_T,Beta_v,Beta_s,Alpha_p,X,Y,Z,t_vec,ke_total,Invariants,[name_file_out '_' num2str(t_out(tt))],b_flag);
            end
        end

%         drawnow
%         plot(P(2,:,2)); hold on


               
     end
catch ME
    ME.message
    disp(['Simulation stopped at time = ' num2str(time) 's'])
    DataOutput(u,v,w,E,rho,mu,kappa,P,T,ke,e,c_p,c_v,sos,Beta_T,Beta_v,Beta_s,Alpha_p,X,Y,Z,t_vec,ke_total,Invariants,[name_file_out '_' num2str(t_out(tt))],b_flag);
end
toc

%% Update last time step data
t = t + 1;
t_vec(t)         = t_vec(t) + time;
ke_total(t)      = Calculate_ke(u, v, w, X, Y, Z);
Invariants       = Calculate_Invariants(dx,dy,dz,rho, u, v, w, E, e, c_v , gamma, T, P, t,Invariants,Ref,Substance,Test,bSolver);

%% Save data
DataOutput(u,v,w,E,rho,mu,kappa,P,T,ke,e,c_p,c_v,sos,Beta_T,Beta_v,Beta_s,Alpha_p,X,Y,Z,t_vec,ke_total,Invariants,[name_file_out '_' num2str(t_out(tt))],b_flag);

% Trim invariants to flag for live processing
[t_vec,ke_total,Invariants] = TrimTimeVariables(t_vec,ke_total,Invariants,b_flag);


%% Print data to file
disp("Simulation completed at t = " + num2str(time) + " s")


end

