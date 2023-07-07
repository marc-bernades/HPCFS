function HP_TransportCoefs = Calculate_HP_TransportCoefs(Substance)

[MW, Tc, pc, p_inf, rhoc, vc_bar, omega, gamma, e_0, c_v, NASA_coefficients, ...
    mu_0, kappa_0, T_0, S_mu, S_kappa, dipole_moment, association_factor] = Substance_Library(Substance);

    % Adimensional dipole moment -- Poling et al. The properties of gases and liquids. McGraw-Hill; 2001.
    HP_TransportCoefs.adimensional_dipole_moment = 131.3*( dipole_moment/sqrt( ( 1.0e6*vc_bar )*Tc ) );
      
    % Viscosity mu -- Poling et al. The properties of gases and liquids. McGraw-Hill; 2001. (9.40; Table 9-6)
    a1_mu  = 6.324;    b1_mu  = 50.412;    c1_mu  = -51.680;   d1_mu  = 1189.0;
    a2_mu  = 1.210e-3; b2_mu  = -1.154e-3; c2_mu  = -6.257e-3; d2_mu  = 0.03728;
    a3_mu  = 5.283;    b3_mu  = 254.209;   c3_mu  = -168.48;   d3_mu  = 3898.0;
    a4_mu  = 6.623;    b4_mu  = 38.096;    c4_mu  = -8.464;    d4_mu  = 31.42;
    a5_mu  = 19.745;   b5_mu  = 7.630;     c5_mu  = -14.354;   d5_mu  = 31.53;
    a6_mu  = -1.900;   b6_mu  = -12.537;   c6_mu  = 4.985;     d6_mu  = -18.15;
    a7_mu  = 24.275;   b7_mu  = 3.450;     c7_mu  = -11.291;   d7_mu  = 69.35;
    a8_mu  = 0.7972;   b8_mu  = 1.117;     c8_mu  = 0.01235;   d8_mu  = -4.117;
    a9_mu  = -0.2382;  b9_mu  = 0.06770;   c9_mu  = -0.8163;   d9_mu  = 4.025;
    a10_mu = 0.06863;  b10_mu = 0.3479;    c10_mu = 0.5926;    d10_mu = -0.727;

    HP_TransportCoefs.E1_mu  = a1_mu  + b1_mu*omega  + c1_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d1_mu*association_factor;
    HP_TransportCoefs.E2_mu  = a2_mu  + b2_mu*omega  + c2_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d2_mu*association_factor;
    HP_TransportCoefs.E3_mu  = a3_mu  + b3_mu*omega  + c3_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d3_mu*association_factor;
    HP_TransportCoefs.E4_mu  = a4_mu  + b4_mu*omega  + c4_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d4_mu*association_factor;
    HP_TransportCoefs.E5_mu  = a5_mu  + b5_mu*omega  + c5_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d5_mu*association_factor;
    HP_TransportCoefs.E6_mu  = a6_mu  + b6_mu*omega  + c6_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d6_mu*association_factor;
    HP_TransportCoefs.E7_mu  = a7_mu  + b7_mu*omega  + c7_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d7_mu*association_factor;
    HP_TransportCoefs.E8_mu  = a8_mu  + b8_mu*omega  + c8_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d8_mu*association_factor;
    HP_TransportCoefs.E9_mu  = a9_mu  + b9_mu*omega  + c9_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0  + d9_mu*association_factor;
    HP_TransportCoefs.E10_mu = a10_mu + b10_mu*omega + c10_mu*HP_TransportCoefs.adimensional_dipole_moment.^4.0 + d10_mu*association_factor;

    % Thermal conductivity k -- Poling et al. The properties of gases and liquids. McGraw-Hill; 2001. (10.23; Table 10-3)
    a1_k = 2.4166*1.00; b1_k = 7.4824*0.100; c1_k = -9.1858*0.10; d1_k = 1.2172*100.0;
    a2_k = -5.0924*0.1; b2_k = -1.5094*1.00; c2_k = -4.9991*10.0; d2_k = 6.9983*10.00;
    a3_k = 6.6107*1.00; b3_k = 5.6207*1.000; c3_k = 6.4760*10.00; d3_k = 2.7039*10.00;
    a4_k = 1.4543*10.0; b4_k = -8.9139*1.00; c4_k = -5.6379*1.00; d4_k = 7.4344*10.00;
    a5_k = 7.9274*0.10; b5_k = 8.2019*0.100; c5_k = -6.9369*0.10; d5_k = 6.3173*1.000;
    a6_k = -5.8634*1.0; b6_k = 1.2801*10.00; c6_k = 9.5893*1.000; d6_k = 6.5529*10.00;
    a7_k = 9.1089*10.0; b7_k = 1.2811*100.0; c7_k = -5.4217*10.0; d7_k = 5.2381*100.0;

    HP_TransportCoefs.B1_k = a1_k + b1_k*omega + c1_k*HP_TransportCoefs.adimensional_dipole_moment.^4.0 + d1_k*association_factor;
    HP_TransportCoefs.B2_k = a2_k + b2_k*omega + c2_k*HP_TransportCoefs.adimensional_dipole_moment.^4.0 + d2_k*association_factor;
    HP_TransportCoefs.B3_k = a3_k + b3_k*omega + c3_k*HP_TransportCoefs.adimensional_dipole_moment.^4.0 + d3_k*association_factor;
    HP_TransportCoefs.B4_k = a4_k + b4_k*omega + c4_k*HP_TransportCoefs.adimensional_dipole_moment.^4.0 + d4_k*association_factor;
    HP_TransportCoefs.B5_k = a5_k + b5_k*omega + c5_k*HP_TransportCoefs.adimensional_dipole_moment.^4.0 + d5_k*association_factor;
    HP_TransportCoefs.B6_k = a6_k + b6_k*omega + c6_k*HP_TransportCoefs.adimensional_dipole_moment.^4.0 + d6_k*association_factor;
    HP_TransportCoefs.B7_k = a7_k + b7_k*omega + c7_k*HP_TransportCoefs.adimensional_dipole_moment.^4.0 + d7_k*association_factor;






end