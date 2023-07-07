%% Time integration check 1D Adv > Rho evolving

clc; clear; close all

N = 20;
h = 1/N;
x = 0:h:(1-h)+h/2;

% y0 = exp(-100*(x-0.5).^2); y0 = y0(:);
y0 = sin(2*pi*x); y0 = y0(:);

dt = 0.001;

[a,b,c] = ButcherTableau(3);

% Derivative matrix for periodic boundaries
v = zeros(1,N); v(2) = 1/2/h; v(end) = -1/2/h;
D = gallery('circul',v);

ys = y0;
for i=1:1/dt

    for i_RK=1:3
        k{i_RK} = -D*ys;
        ys = TimeIntegration_RKGeneral(y0,dt,k,i_RK,3,a,b,c);
    end
    y0 = ys;
%     plot(y0)
%     drawnow
%     pause(0.01)

end
