function [mu,kappa] = Calculate_muandkappa(HP_TransportCoefs,Substance, T, rho)

% Convert HP_Transport Coefficients into variables
f = fieldnames(HP_TransportCoefs);
for index = 1:length(f)
    eval([f{index} ' = HP_TransportCoefs. ' f{index} ';']);
end

% Substance
[MW, Tc, pc, p_inf, rhoc, vc_bar, omega, gamma, e_0, c_v, NASA_coefficients, ...
    mu_0, kappa_0, T_0, S_mu, S_kappa, dipole_moment, association_factor] = Substance_Library(Substance);

acentric_factor = omega;

%% Dynamic viscosity
% Chung et al 1984 and 1988
v          = MW./rho;
Y          = vc_bar./(6*v);
T_ast      = 1.2593*( T/Tc );
Omega      = 1.16145*(T_ast.^-0.14874  + 0.52487*exp(-0.77320*T_ast) + 2.16178*exp(-2.43787*T_ast));
G1         = ( 1.0 - 0.5*Y )./(1.0 - Y).^3.0;
Fc         = 1.0 - 0.2756*acentric_factor + 0.059035*( adimensional_dipole_moment.^4.0 ) + association_factor;
G2_mu      = ( E1_mu*( 1.0 - exp( -E4_mu*Y ) )./Y + E2_mu*G1.*exp( E5_mu*Y ) + E3_mu*G1 )./( E1_mu*E4_mu + E2_mu + E3_mu );
mu_ast_ast = E7_mu*( Y.^2.0 ).*G2_mu.*exp( E8_mu + E9_mu./T_ast + E10_mu.*( T_ast.^-2.0 ) );
mu_ast     = ( sqrt( T_ast ).*Fc./Omega ).*( 1.0./G2_mu + E6_mu.*Y ) + mu_ast_ast;

mu         = 1.0e-7*mu_ast.*( ( 36.344*sqrt( ( 1.0E3*MW ).*Tc ) )./( (1.0e6*vc_bar).^(2.0/3.0) ) );

%% Thermal conductivity viscosity
% PengRobinson
[a,b,R,dadT,d2adT2,NASA_coefficients] = PengRobinson(T, Substance);

% Calculate ideal Cp depending on temperature
c_p_ideal = zeros(size(T)); % Preallocate
for i = 1:length(T(:,1,1))
    for j = 1:length(T(1,:,1))
        for k = 1:length(T(1,1,:))
            if T(i,j,k) >= 200 && T(i,j,k) < 1000
                c_p_ideal(i,j,k) = R.*(NASA_coefficients(8) + NASA_coefficients(9)*T(i,j,k) + NASA_coefficients(10)*T(i,j,k).^2 + NASA_coefficients(11)*T(i,j,k).^3 + NASA_coefficients(12)*T(i,j,k).^4);
            elseif T(i,j,k) >= 1000 && T(i,j,k) < 6000
                c_p_ideal(i,j,k) = R.*(NASA_coefficients(1) + NASA_coefficients(2)*T(i,j,k) + NASA_coefficients(3)*T(i,j,k).^2 + NASA_coefficients(4)*T(i,j,k).^3 + NASA_coefficients(5)*T(i,j,k).^4);
            elseif T(i,j,k) < 200
                % Assume constant temperature below 200K
                c_p_ideal(i,j,k) = R.*(NASA_coefficients(8) + NASA_coefficients(9)*200 + NASA_coefficients(10)*200.^2 + NASA_coefficients(11)*200.^3 + NASA_coefficients(12)*200.^4);
            else
                break
            end
        end
    end
end
% Molar enthalpy
std_bar_c_p = c_p_ideal*MW;
R_universal = R*MW;

% Chung et al 1984 and 1988
mu_0_k  = 40.785e-7*Fc.*sqrt( 1.0e3*MW.*T )./( ( (1.0e6*vc_bar).^(2.0/3.0))*Omega );
alpha_k = ( std_bar_c_p/R_universal - 1.0 ) - 1.5;
beta_k  = 0.7862 - 0.7109*acentric_factor + 1.3168*( acentric_factor^2.0 );
gamma_k = 2.0 + 10.5*( (T/Tc).^2.0 );
Psi_k   = 1.0 + alpha_k.*( ( 0.215 + 0.28288*alpha_k - 1.061*beta_k + 0.26665*gamma_k )./( 0.6366 + beta_k*gamma_k + 1.061*alpha_k*beta_k) );
q_k     = 0.003586*( sqrt( Tc/MW )./( ( 1.0e6*vc_bar ).^(2.0/3.0) ) );
G3_k    = ( ( ( B1_k./Y ).*( 1.0 - exp( (-1.0)*B4_k.*Y ) ) ) + ( B2_k.*G1.*exp( B5_k.*Y ) ) + ( B3_k.*G1 ) )./( B1_k.*B4_k + B2_k + B3_k );

% Calculate thermal conductivity
kappa   = ( 31.2*mu_0_k.*Psi_k./MW ).*( 1.0./G3_k + B6_k*Y ) + q_k*B7_k.*(Y.^2.0).*sqrt( T./Tc).*G3_k;




end
