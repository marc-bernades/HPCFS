function [ a,b,R,dadT,d2adT2,NASA_coefficients ] = PengRobinson( T, Substance )

% getThermo Compute necessary thermodyamic parameters for the cubic
% equation of state (N2 property is assumed).

% N2 substance library
[MW, Tc, pc, p_inf, rhoc, vc_bar, omega, gamma, e_0, c_v, NASA_coefficients, ...
    mu_0, kappa_0, T_0, S_mu, S_kappa, dipole_moment, association_factor] = Substance_Library(Substance);



% Peng Robsinson coefficients
if omega > 0.49 % Accentric factor
    c      = 0.379642 + 1.48503*omega - 0.164423*omega^2 + 0.016666*omega^3;
else
    c      = 0.37464 + 1.54226*omega - 0.26992*omega^2;
end

Ru     = 8.314;   % R universal
R      = Ru/MW;   % R specific
a      = 0.457236*(R*Tc)^2/pc*(1+c*(1-sqrt(T/Tc))).^2;
b      = 0.077796*R*Tc./pc;
G      = c*sqrt(T/Tc)./(1+c*(1-sqrt(T/Tc)));
dadT   = -1./T.*a.*G;
d2adT2 = 0.457236*R^2./T/2*c*(1+c)*Tc/pc.*sqrt(Tc./T);

end

