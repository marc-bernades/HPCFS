function [ rho ] = Calculate_Rho_from_TandP( bSolver, T,p, Fluid, Substance )


if strcmp(bSolver,'Real')
    % PengRobinson
    [a,b,R,dadT,d2adT2,NASA_coefficients] = PengRobinson(T, Substance);

    % Auxiliar parameters
    a0 = p*b^3 + b^2*R*T - a*b;
    a1 = -3*p*b^2 - 2*b*R*T + a;
    a2 = p*b-R*T;
    a3 = p;

    % Only one real roots
    l = a2./a3;
    m = a1./a3;
    n = a0./a3;
    A = 1/3*(3*m-l.^2);
    B = 1/27*(2*l.^3 - 9*l.*m + 27*n);
    D = A.^3/27 + B.^2/4;
    M = sign(-B/2+sqrt(D)).*(abs(-B/2+sqrt(D))).^(1/3);
    N = sign(-B/2-sqrt(D)).*(abs(-B/2-sqrt(D))).^(1/3);

    v = M + N - l/3;
    rho = 1./v;

%     rho =  p./(Fluid.R_specific*T);

else % IDEAL GAS
    rho   = p./(Fluid.R_specific*T);

end

end

