function flux  = WENO5_solver(Fp,Fn)
%%  WENO5 setup parameters
C_p     = [2/6 -7/6 11/6 0    0;
            0  -1/6 5/6 1/3   0;
            0    0  1/3 5/6 -1/6];  % Positive matrix for polynomial
C_n     = flip(flip(C_p),2);        % Negative matrix for -ve part
d       = [1/10 6/10 3/10];         % Gamma vector
Epsilon = 1E-6;


%% Positive Flux > Extrapolation $v_{i+1/2}^{-}$ == $f_{i+1/2}^{+}$
vmm = Fn(:,1);
vm  = Fn(:,2);
v0  = Fn(:,3);
vp  = Fn(:,4);
vpp = Fn(:,5);

% Smooth indicators Beta factors
Beta_0 = 13/12*(vmm - 2*vm + v0).^2  + 1/4*(vmm  - 4*vm + 3*v0).^2;
Beta_1 = 13/12*(vm  - 2*v0 + vp).^2  + 1/4*(vm   -   vp).^2;
Beta_2 = 13/12*(v0  - 2*vp + vpp).^2 + 1/4*(3*v0 - 4*vp + vpp).^2;

% Alpha weights
alpha_0 = d(1)./(Epsilon + Beta_0).^2;
alpha_1 = d(2)./(Epsilon + Beta_1).^2;
alpha_2 = d(3)./(Epsilon + Beta_2).^2;

% ENO Stencil weights
w0 = alpha_0./sum(alpha_0 + alpha_1 + alpha_2);
w1 = alpha_1./sum(alpha_0 + alpha_1 + alpha_2);
w2 = alpha_2./sum(alpha_0 + alpha_1 + alpha_2);


% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
flux_p = w0.*(C_p(1,1)*vmm + C_p(1,2)*vm + C_p(1,3)*v0 + C_p(1,4)*vp + C_p(1,5)*vpp) ...
    + w1.*(C_p(2,1)*vmm + C_p(2,2)*vm + C_p(2,3)*v0 + C_p(2,4)*vp + C_p(2,5)*vpp) ...
    + w2.*(C_p(3,1)*vmm + C_p(3,2)*vm + C_p(3,3)*v0 + C_p(3,4)*vp + C_p(3,5)*vpp);


%% Nevgative flux > Extrapolation $u_{i+1/2}^{+}$ == $f_{i+1/2}^{-}$
umm = Fp(:,1);
um  = Fp(:,2);
u0  = Fp(:,3);
up  = Fp(:,4);
upp = Fp(:,5);

% Smooth indicators Beta factors
Beta_0 = 13/12*(umm - 2*um + u0).^2  + 1/4*(umm  - 4*um + 3*u0).^2;
Beta_1 = 13/12*(um  - 2*u0 + up).^2  + 1/4*(um   -   up).^2;
Beta_2 = 13/12*(u0  - 2*up + upp).^2 + 1/4*(3*u0 - 4*up + upp).^2;

% Alpha weights
alpha_0 = d(1)./(Epsilon + Beta_0).^2;
alpha_1 = d(2)./(Epsilon + Beta_1).^2;
alpha_2 = d(3)./(Epsilon + Beta_2).^2;

% ENO Stencil weights
w0 = alpha_0./(alpha_0 + alpha_1 + alpha_2);
w1 = alpha_1./(alpha_0 + alpha_1 + alpha_2);
w2 = alpha_2./(alpha_0 + alpha_1 + alpha_2);


% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
flux_n = w0.*(C_n(1,1)*umm + C_n(1,2)*um + C_n(1,3)*u0 + C_n(1,4)*up + C_n(1,5)*upp) ...
    + w1.*(C_n(2,1)*umm + C_n(2,2)*um + C_n(2,3)*u0 + C_n(2,4)*up + C_n(2,5)*upp) ...
    + w2.*(C_n(3,1)*umm + C_n(3,2)*um + C_n(3,3)*u0 + C_n(3,4)*up + C_n(3,5)*upp);


%% Final flux
flux = flux_p + flux_n;




end
