function flux  = WENO2_solver(Fp,Fn)
%%  WENO5 setup parameters

C_p     = [-1/2 3/2  0;
             0  1/2 1/2];  % Positive matrix for polynomial
C_n     = flip(flip(C_p),2);        % Negative matrix for -ve part
d       = [1/3 2/3];         % Gamma vector

Epsilon = 1E-6;


%% Positive Flux > Extrapolation $v_{i+1/2}^{-}$ == $f_{i+1/2}^{+}$
vm  = Fp(:,1);
v0  = Fp(:,2);
vp  = Fp(:,3);

% Smooth indicators Beta factors
Beta_0 = (vm - v0).^2;
Beta_1 = (v0 - vp).^2;

% Alpha weights
alpha_0 = d(1)./(Epsilon + Beta_0).^2;
alpha_1 = d(2)./(Epsilon + Beta_1).^2;

% ENO Stencil weights
w0 = alpha_0./sum(alpha_0 + alpha_1);
w1 = alpha_1./sum(alpha_0 + alpha_1);

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
flux_p = w0.*(C_p(1,1)*vm + C_p(1,2)*v0 + C_p(1,3)*vp) ...
    + w1.*(C_p(2,1)*vm + C_p(2,2)*v0 + C_p(2,3)*vp);

%% Nevative flux > Extrapolation $u_{i+1/2}^{+}$ == $f_{i+1/2}^{-}$
um = Fn(:,1);
u0  = Fn(:,2);
up  = Fn(:,3);

% Smooth indicators Beta factors
Beta_0 = (um - u0).^2;
Beta_1 = (u0 - up).^2;

% Alpha weights
alpha_0 = d(1)./(Epsilon + Beta_0).^2;
alpha_1 = d(2)./(Epsilon + Beta_1).^2;

% ENO Stencil weights
w0 = alpha_0./(alpha_0 + alpha_1);
w1 = alpha_1./(alpha_0 + alpha_1);


% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
flux_n = w0.*(C_n(1,1)*um + C_n(1,2)*u0 + C_n(1,3)*up) ...
    + w1.*(C_n(2,1)*um + C_n(2,2)*u0 + C_n(2,3)*up);


%% Final flux
flux = flux_p + flux_n;




end
