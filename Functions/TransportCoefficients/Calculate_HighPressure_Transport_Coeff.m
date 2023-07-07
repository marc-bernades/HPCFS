function [mu,kappa] = Calculate_HighPressure_Transport_Coeff(bSolver, mu, kappa, Fluid, Substance, T, rho, HP_model)

if strcmp(bSolver,'Real')
    % Substance
    [MW, Tc, pc, p_inf, rhoc, vc_bar, omega, gamma, e_0, c_v, NASA_coefficients, ...
        mu_0, kappa_0, T_0, S_mu, S_kappa, dipole_moment, association_factor] = Substance_Library(Substance);

    switch HP_model
        case 'Constant'
            %% Constant model low pressure
            if isfield(Fluid,'mu_0') && isfield(Fluid,'kappa_0')
                % User definition
                mu    = zeros(size(mu)) + Fluid.mu_0;
                kappa = zeros(size(kappa)) + Fluid.kappa_0;

            else
                %From substance library
                mu    = zeros(size(mu))    + mu_0;
                kappa = zeros(size(kappa)) + kappa_0;
            end

        case 'LowPressure'

            %% Low pressure transport coefficients Sutherland's law of viscosity
            mu    = (mu_0.*(T/T_0).^1.5).*(T_0 + S_mu)./(T + S_mu);
            kappa = (kappa_0.*(T/T_0).^1.5).*(T_0 + S_kappa)./(T + S_kappa);

        case 'HighPressure'

            %% High pressure model
            HP_TransportCoefs = Calculate_HP_TransportCoefs(Substance);
            [mu,kappa] = Calculate_muandkappa(HP_TransportCoefs,Substance, T, rho);

    end


else

    switch HP_model
        case 'Constant'
            %Ideal gas > From substance library
            mu    = zeros(size(mu)) + Fluid.mu_0;
            kappa = zeros(size(kappa)) + Fluid.kappa_0;

        case 'Power'

            % Reference Sutherland law for IG
            T_0     = 273;       % K
            mu_0    = 1.37*1E-5; % Pa s
            kappa_0 = 0.0146;    % W / (m K)
            n1      = 0.79;
            n2      = 1.30;

            %% Low pressure transport coefficients Power law of viscosity
            mu    = (mu_0.*(T/T_0).^n1);
            kappa = (kappa_0.*(T/T_0).^n2);



        case 'LowPressure'

            % Reference Sutherland law for IG
            T_0     = 273;       % K
            mu_0    = 1.37*1E-5; % Pa s
            kappa_0 = 0.0146;    % W / (m K)
            S_mu    = 222;       % K
            S_kappa = 1800;      % K

            %% Low pressure transport coefficients Sutherland's law of viscosity
            mu    = (mu_0.*(T/T_0).^1.5).*(T_0 + S_mu)./(T + S_mu);
            kappa = (kappa_0.*(T/T_0).^1.5).*(T_0 + S_kappa)./(T + S_kappa);

    end


end

end