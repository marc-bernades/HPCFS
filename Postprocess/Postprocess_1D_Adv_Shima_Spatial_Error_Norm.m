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

t_target          = 0.5*Fluid.t_final;       % Target time to show invariants

Scheme_label = {'KGP-Et','KGP-Et-CR','KGP-et-CR','KGP-Pt'};

% KGP-Et
Scheme{1,1}       = readtable('output_data_1D_Adv_Shima_20x1x1_KGP-Et_0.01_Time.csv');
Scheme{1,2}       = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-Et_0.01_Time.csv');
Scheme{1,3}       = readtable('output_data_1D_Adv_Shima_82x1x1_KGP-Et_0.01_Time.csv');
Scheme{1,4}       = readtable('output_data_1D_Adv_Shima_164x1x1_KGP-Et_0.01_Time.csv');

% KGP-Et-CR
Scheme{2,1}       = readtable('output_data_1D_Adv_Shima_20x1x1_KGP-Et-CR_0.01_Time.csv');
Scheme{2,2}       = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-Et-CR_0.01_Time.csv');
Scheme{2,3}       = readtable('output_data_1D_Adv_Shima_82x1x1_KGP-Et-CR_0.01_Time.csv');
Scheme{2,4}       = readtable('output_data_1D_Adv_Shima_164x1x1_KGP-Et-CR_0.01_Time.csv');

% KGP-et-CR
Scheme{3,1}       = readtable('output_data_1D_Adv_Shima_20x1x1_KGP-et-CR_0.01_Time.csv');
Scheme{3,2}       = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-et-CR_0.01_Time.csv');
Scheme{3,3}       = readtable('output_data_1D_Adv_Shima_82x1x1_KGP-et-CR_0.01_Time.csv');
Scheme{3,4}       = readtable('output_data_1D_Adv_Shima_164x1x1_KGP-et-CR_0.01_Time.csv');

% KGP-Pt
Scheme{4,1}       = readtable('output_data_1D_Adv_Shima_20x1x1_KGP-Pt_0.01_Time.csv');
Scheme{4,2}       = readtable('output_data_1D_Adv_Shima_41x1x1_KGP-Pt_0.01_Time.csv');
Scheme{4,3}       = readtable('output_data_1D_Adv_Shima_82x1x1_KGP-Pt_0.01_Time.csv');
Scheme{4,4}       = readtable('output_data_1D_Adv_Shima_164x1x1_KGP-Pt_0.01_Time.csv');

%% Iterate schemes for spatial error plot on PEP / TEC
for i = 1:length(Scheme)

    for j = 1:length(N)
        % Index only unique Y and Z
        try
            [val,Index]=min(abs(Scheme{i,j}.t_vec-t_target));
%             Index=Scheme{i,j}.t_vec(idx);
%             [rho_ref, P_ref, T_ref, e_ref, u_ref, v_ref, w_ref, ke_ref, E_ref ] = Reference_Solution_1D_Adv_Shima(bSolver,X,Fluid,Substance);

            % PEP error
            %         [P_L1_norm_error(i,j), ~]    = Calculate_Error_L1(Scheme{i,j}.P(Index),[],P_exact,[],h_x(j),[],[],Test);
            P_L2_norm_error(i,j)    = abs(Scheme{i,j}.P_norm(Index))*Fluid.P_0;

            % TEC error
            %         [E_L1_norm_error(i,j), ~]    = Calculate_Error_L1(Scheme{i,j}.E(Index),[],E_ref,[],h_x(j),[],[],Test);
            %         [E_L2_norm_error(i,j), ~]    = Calculate_Error_L2(Scheme{i,j}.E(Index),[],E_ref,[],h_x(j),[],[],Test);
%             rhoE_ref = sum(sum(sum(rho_ref.*E_ref)));
            E_L2_norm_error(i,j) = Scheme{i,j}.rhoE_norm(Index);


        catch
            P_L2_norm_error(i,j) = nan;
            E_L2_norm_error(i,j) = nan;
        end
    end


    %     plot_SpatialError(P_L1_norm_error,[], P_L2_norm_error,[],L_x,[],[],N,Test)
    %     plot_SpatialError(E_L1_norm_error,[], E_L2_norm_error,[],L_x,[],[],N,Test)


end

% Colors KGP-Et, KGP-Et-CR, KGP-et-CR, KGP-Pt
cSchemes = {[0.3010 0.7450 0.9330],[0.3010 0.7450 0.9330],[0, 0, 1],[0 0.4470 0.7410]};
styleSchemes = {'--',':','-.','-'};

% Plots
% PEP Plot
figure
% subplot(1,2,1)
offset_Oh2_PEP = 10^3;
for i = 1:length(Scheme)-1
    loglog(L_x./N,abs(P_L2_norm_error(i,:)),'LineWidth',2, 'LineStyle',styleSchemes{i},'Marker','o','MarkerSize',5,'color',cSchemes{i});
    hold on; grid on;
end
loglog(L_x./N,(offset_Oh2_PEP*L_x./N.^2),'k--')
% xtickformat('10^-^%d')
legend([Scheme_label(1:end-1),'${h}^2$'],'interpreter','latex','fontsize',12,'Location','northwest')
ylabel('PEP L2 Norm')
legend('Location','northwest','box','off')
% pbaspect([1 1.1 1])
% set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{h} ({L}/{N})}$','interpreter','latex','fontsize',16)
ylabel('$\varepsilon_{PEP}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\1D_Adv_Shima\Real\41x1x1\1D_Shima_Real_Error_PEP','epsc')



% TEC Plot
figure
offset_Oh2_TEC = 10^-7;
for i = 1:length(Scheme)
    loglog(L_x./N,abs(E_L2_norm_error(i,:)),'LineWidth',2, 'LineStyle',styleSchemes{i},'Marker','o','MarkerSize',5,'color',cSchemes{i});
    hold on; grid on;
end
loglog(L_x./N,(offset_Oh2_TEC*L_x./N.^2),'k--')
% xtickformat('10^-^%d')
legend([Scheme_label,'${h}^2$'],'interpreter','latex','fontsize',12,'Location','northwest')
ylabel('PEP L2 Norm')
legend('Location','northwest','box','off')
% pbaspect([1 1.1 1])
% set(gca,'linewidth',1.5)
set(gca,'fontsize',12)
xlabel('${{h} ({L}/{N})}$','interpreter','latex','fontsize',16)
ylabel('$\varepsilon_{TEC}$','interpreter','latex','fontsize',16)
saveas(gca,'Results\1D_Adv_Shima\Real\41x1x1\1D_Shima_Real_Error_TEC','epsc')