function [u_pert,v_pert] = Calculate_2D_MixingLayer_Perturbation(u,dy,X,Y,L_y)

% Free stream top and bottom speed
U_inf_top = 0.5*(u(end,2,2) + u(end-1,2,2));
U_inf_bot = 0.5*(u(1,2,2)   + u(2,2,2));
U_inf_ref = U_inf_top - U_inf_bot;            % Reference speed removing free-slip velocity
% Momentum thickness (ref https://www.sciencedirect.com/topics/engineering/momentum-thickness)
theta = 0;
% Integrate entire domain in y-direction
for ii = 1:length(u(:,2,2))
    theta = theta + dy(ii,2,2)*((u(ii,2,2) - U_inf_bot)/U_inf_ref.*(1 - (u(ii,2,2) - U_inf_bot)/U_inf_ref));
end

% Wave number (ref Sharman JFM 2019)
k      = 0.22/theta; % Wave number most unstable mode
k      = 6; % To mimic similar of pi*k than 2*k (ref) to set it within domain
lambda = 2*pi/(k*pi);     % Wave number should be lower than the domain


% Perturbation amplitude
A   = 0.1;% maximum velocity perturbation amplitude of approximately 7.5% of the free-stream flow-speed difference across the shear layer
A_g = 1/10;
a_y = U_inf_ref.*A*(exp(-((Y - L_y/2).^2/theta)) + exp(-(Y.^2/theta)) + exp(-((Y + L_y/2).^2/theta)));

% Perturbation fields
u_pert = a_y.*(sin(k*pi*X) + A_g*rand(size(X))/2);
v_pert = a_y.*(sin(k*pi*X) + A_g*rand(size(X))/2);

end