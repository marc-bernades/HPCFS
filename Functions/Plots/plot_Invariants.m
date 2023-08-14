function plot_Invariants(t_vec, rho_norm, rhou_bar, rhoE_norm, rhoe_norm, T_norm, rhoui2_norm, ...
    rhoE2_norm, rhoe2_norm, P_norm, varargin)

if isempty(varargin)
    % Only Solver results
    figure()
    subplot(3,3,1)
    plot(t_vec, rho_norm)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho} \rangle}$','interpreter','latex')

    subplot(3,3,2)
    plot(t_vec, rhou_bar)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{{\rho}{u}}$','interpreter','latex')

    subplot(3,3,3)
    plot(t_vec, rhoE_norm)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{E} \rangle}$','interpreter','latex')

    subplot(3,3,4)
    plot(t_vec, rhoe_norm)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{e} \rangle}$','interpreter','latex')

    subplot(3,3,5)
    plot(t_vec, T_norm)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{\langle {T} \rangle}$','interpreter','latex')

    subplot(3,3,6)
    plot(t_vec, rhoui2_norm)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{{u_i}^2} \rangle}$','interpreter','latex')

    subplot(3,3,7)
    plot(t_vec, rhoE2_norm)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{{E}^2} \rangle}$','interpreter','latex')

    subplot(3,3,8)
    plot(t_vec, rhoe2_norm)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{{e}^2} \rangle}$','interpreter','latex')

    subplot(3,3,9)
    plot(t_vec, P_norm)
    xlim([0 max(t_vec)])
    xlabel('time [s]')
    ylabel('$\overline{\langle {P} \rangle}$','interpreter','latex')

else
    %JCP_ref = {t, rho_F, rhou_F, rhoE_F, rhoe_F, rhos_F, rhouq_F, rhoEq_F, rhoeq_F, rhosq_F};
    JCP_ref = {'t', 'rho_F', 'rhou_F', 'rhoE_F', 'rhoe_F', 'rhos_F', 'rhouq_F', 'rhoEq_F', 'rhoeq_F', 'rhosq_F'};

    for ii = 1:length(JCP_ref)
        % Create variables in structure
        JCP.([JCP_ref{ii}]) = varargin{ii};
        % Set them in workspace
        eval([JCP_ref{ii} '=JCP.' JCP_ref{ii} ]);
    end

    % Plot results
    figure()
    subplot(3,3,1)
    plot(t_vec, rho_norm)
    hold on
    plot(t,rho_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho} \rangle}$','interpreter','latex')

    subplot(3,3,2)
    plot(t_vec, rhou_bar)
    hold on
    plot(t,rhou_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{{\rho}{u}}$','interpreter','latex')

    subplot(3,3,3)
    plot(t_vec, rhoE_norm)
    hold on
    plot(t,rhoE_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{E} \rangle}$','interpreter','latex')

    subplot(3,3,4)
    plot(t_vec, rhoe_norm)
    hold on
    plot(t,rhoe_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{e} \rangle}$','interpreter','latex')

    subplot(3,3,5)
    plot(t_vec, rhos_norm)
    hold on
    plot(t,rhos_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{s} \rangle}$','interpreter','latex')

    subplot(3,3,6)
    plot(t_vec, rhoui2_norm)
    hold on
    plot(t,rhouq_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{{u_i}^2} \rangle}$','interpreter','latex')

    subplot(3,3,7)
    plot(t_vec, rhoE2_norm)
    hold on
    plot(t,rhoEq_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{{E}^2} \rangle}$','interpreter','latex')

    subplot(3,3,8)
    plot(t_vec, rhoe2_norm)
    hold on
    plot(t,rhoeq_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{{e}^2} \rangle}$','interpreter','latex')

    subplot(3,3,9)
    plot(t_vec, rhos2_norm)
    hold on
    plot(t,rhosq_F)
    xlim([0 max(t_vec)])
    legend('Solver','JCP ref')
    xlabel('time [s]')
    ylabel('$\overline{\langle {\rho}{{s}^2} \rangle}$','interpreter','latex')

end



end