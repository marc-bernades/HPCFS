function X = TimeIntegration_HighStabilityRK3(X,X_0,dt,k,i_RK)
c   = [0, 1/2, 1/2, 1];
a   = [[0, 0, 0, 0]
       [1/2, 0, 0, 0]
       [0, 1/2, 0, 0]
       [0, 0, 1, 0]];
b   = [1/6, 1/3, 1/3, 1/6];

if i_RK == 1
    X = X_0 + dt*k{i_RK};
elseif i_RK == 2
    X = X_0 + dt/4*(k{i_RK-1} + k{i_RK});
elseif i_RK == 3
    X = X_0 + dt/6*(k{i_RK-2} + k{i_RK - 1} + 4*k{i_RK});

end

end

