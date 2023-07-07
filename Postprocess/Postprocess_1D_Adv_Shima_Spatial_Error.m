%% Postprocess PEP error Adv
bSolver           = 'Real';
Substance         = 'N2';                      % Fluid selected substance
HP_model          = 'Constant';                % Transport coefficients model: 'Constant', 'LowPressure', 'HighPressure'

Test     = '1D_Adv';
P_exact  = 5E6;
L_x      = 1;
N        = [20 41 82 164];
h_x      = L_x./N;

Fluid.U_0         = 1;                      % Reference velocity [m/s]
Fluid.P_0         = 5E6;                    % Rereference Pressure [Pa]
Fluid.rho_min     = 56.9;                   % Min density non-linear
Fluid.rho_max     = 793.1;                  % Max density non-linear
Fluid.Case        = 'Smooth';               % Non-linear case "Smooth" as per Shima et al.
Fluid.mu_0        = 0;
Fluid.kappa_0     = 0;
Fluid.t_final     = 0.01;

Scheme_label = {'KGP-Et','KGP-Et-CR','KGP-et-CR','KGP-Pt'};

% KGP-Et
Scheme{1,1}       = readtable('output_data_1D_Adv_Shima_20x1x1_KGP-Et_0.01.csv');
Scheme{1,2}       = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-Et_0.01.csv');
Scheme{1,3}       = readtable('output_data_1D_Adv_Shima_82x1x1_KGP-Et_0.01.csv');
Scheme{1,4}       = readtable('output_data_1D_Adv_Shima_164x1x1_KGP-Et_0.01.csv');

% KGP-Et-CR
Scheme{2,1}       = readtable('output_data_1D_Adv_Shima_20x1x1_KGP-Et-CR_0.01.csv');
Scheme{2,2}       = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-Et-CR_0.01.csv');
Scheme{2,3}       = readtable('output_data_1D_Adv_Shima_82x1x1_KGP-Et-CR_0.01.csv');
Scheme{2,4}       = readtable('output_data_1D_Adv_Shima_164x1x1_KGP-Et-CR_0.01.csv');

% KGP-et-CR
Scheme{3,1}       = readtable('output_data_1D_Adv_Shima_20x1x1_KGP-et-CR_0.01.csv');
Scheme{3,2}       = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-et-CR_0.01.csv');
Scheme{3,3}       = readtable('output_data_1D_Adv_Shima_82x1x1_KGP-et-CR_0.01.csv');
Scheme{3,4}       = readtable('output_data_1D_Adv_Shima_164x1x1_KGP-et-CR_0.01.csv');

% KGP-Pt
Scheme{4,1}       = readtable('output_data_1D_Adv_Shima_20x1x1_KGP-Pt_0.01.csv');
Scheme{4,2}       = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-Pt_0.01.csv');
Scheme{4,3}       = readtable('output_data_1D_Adv_Shima_82x1x1_KGP-Pt_0.01.csv');
Scheme{4,4}       = readtable('output_data_1D_Adv_Shima_164x1x1_KGP-Pt_0.01.csv');

%% Iterate schemes for spatial error plot on PEP / TEC
for i = 1:length(Scheme)

    for j = 1:length(N)
        % Index only unique Y and Z
        try
            Index_Z = unique(Scheme{i,j}.Z); Index_Z = find(Scheme{i,j}.Z == Index_Z(2));
            Index   = Index_Z(5:3:end-3);
            X = Scheme{i,j}.X(Index);
            [rho_ref, P_ref, T_ref, e_ref, u_ref, v_ref, w_ref, ke_ref, E_ref ] = Reference_Solution_1D_Adv_Shima(bSolver,X,Fluid,Substance);

            % PEP error
            %         [P_L1_norm_error(i,j), ~]    = Calculate_Error_L1(Scheme{i,j}.P(Index),[],P_exact,[],h_x(j),[],[],Test);
            [P_L2_norm_error(i,j), ~]    = Calculate_Error_L2(Scheme{i,j}.P(Index),[],P_exact,[],h_x(j),[],[],Test);

            % TEC error
            %         [E_L1_norm_error(i,j), ~]    = Calculate_Error_L1(Scheme{i,j}.E(Index),[],E_ref,[],h_x(j),[],[],Test);
            %         [E_L2_norm_error(i,j), ~]    = Calculate_Error_L2(Scheme{i,j}.E(Index),[],E_ref,[],h_x(j),[],[],Test);
            rhoE_ref = sum(sum(sum(rho_ref.*E_ref)));
            E_L2_norm_error(i,j) = abs((sum(sum(sum(Scheme{i,j}.rho(Index).*Scheme{i,j}.E(Index)))) - rhoE_ref)./rhoE_ref);


        catch
            P_L2_norm_error(i,j) = nan;
            E_L2_norm_error(i,j) = nan;
        end
    end


    %     plot_SpatialError(P_L1_norm_error,[], P_L2_norm_error,[],L_x,[],[],N,Test)
    %     plot_SpatialError(E_L1_norm_error,[], E_L2_norm_error,[],L_x,[],[],N,Test)


end


% Plots
% PEP Plot
figure
subplot(1,2,1)
offset_Oh2_PEP = 10^5.5;
for i = 1:length(Scheme)-1
    loglog(L_x./N,P_L2_norm_error(i,:),'o-','MarkerSize',5); hold on; grid on;
end
loglog(L_x./N,(offset_Oh2_PEP*L_x./N.^2),'k--')
% xtickformat('10^-^%d')
legend([Scheme_label(1:end-1),'h^2'],'Location','northwest')
ylabel('PEP L2 Norm')
xlabel('h (L/N)')

% TEC Plot
subplot(1,2,2)
offset_Oh2_TEC = 10^-8;
for i = 1:length(Scheme)
    loglog(L_x./N,E_L2_norm_error(i,:),'o-','MarkerSize',5); hold on; grid on;
end
loglog(L_x./N,(offset_Oh2_TEC*L_x./N.^2),'k--')
% xtickformat('10^-^%d')
legend([Scheme_label,'h^2'],'Location','northwest')
ylabel('$\overline{\langle {\rho}{E} \rangle}$','interpreter','latex')
xlabel('h (L/N)')

%% Figure Comparison
figure; 
% Select scheme 1: 
i = 1; j = 2;
plot(Scheme{i,j}.X(Index), Scheme{i,j}.P(Index));
hold on
% Select scheme 2:
i = 2; j = 2;
plot(Scheme{i,j}.X(Index), Scheme{i,j}.P(Index));
legend('KGP-Et','KGP-Pt')


%% Exact energy vs Discretized results on the KGP-Pt
figure
for i = 1:length(N)
    % Index only unique Y and Z
    Index_Z = unique(Scheme{i,j}.Z); Index_Z = find(Scheme{i,j}.Z == Index_Z(2));
    Index   = Index_Z(5:3:end-3);
    X = Scheme{i,j}.X(Index);
    [rho_ref, P_ref, T_ref, e_ref, u_ref, v_ref, w_ref, ke_ref, E_ref ] = Reference_Solution_1D_Adv_Shima(bSolver,X,Fluid,Substance);
    
    
    plot(Scheme{i,j}.X(Index), Scheme{i,j}.E(Index)); hold on

    % Exact solution at Nmax
    if i == 3
        plot(X,E_ref,'k--'); hold on
        legend('N = 20','N = 41', 'N = 82', 'Exact')
    end
end
