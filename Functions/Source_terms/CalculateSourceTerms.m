function [f_rhou, f_rhov, f_rhow, f_rhoE, f_rhouvw] = CalculateSourceTerms(f_rhou, f_rhov, f_rhow, f_rhoE,u, v, w, X, Y, Z, K, Fluid)



%% Check controller
if K.b_active == 1

    % Case for transcritical channel flow
    switch K.case
        case '2D_ChannelFlow'
        % Calculate delta_y
        delta_y = Y(2,1,1) - Y(1,1,1);

        % Calculate avg velocities 
        num_points_bw     = max(size(u(:,1,1)));
        avg_u_inner_bw    = sum(sum(sum(u(2,:,:))))/num_points_bw;
        avg_u_boundary_bw = sum(sum(sum(u(1,:,:))))/num_points_bw;

        % Calculate tau_bw_numerical
        tau_bw_numerical = K.mu_bw*( avg_u_inner_bw - avg_u_boundary_bw )/delta_y;

        % Update controller variables
        Controller_error   = (K.tau_bw_target - tau_bw_numerical )/Fluid.delta;
        Controller_output  = K.kp*Controller_error; % Estimated uniform body force to drive the flow


        f_rhou(2:end-1,2:end-1,2:end-1) = Controller_output;
        f_rhov(2:end-1,2:end-1,2:end-1) = 0;
        f_rhow(2:end-1,2:end-1,2:end-1) = 0;
        f_rhoE(2:end-1,2:end-1,2:end-1) = 0;

    end
else
    %% Alternatively Source terms 0
    f_rhou(2:end-1,2:end-1,2:end-1) = 0;
    f_rhov(2:end-1,2:end-1,2:end-1) = 0;
    f_rhow(2:end-1,2:end-1,2:end-1) = 0;
    f_rhoE(2:end-1,2:end-1,2:end-1) = 0;
end


% Calculate momentum sources
f_rhouvw = f_rhou.*u + f_rhov.*v + f_rhow.*w;

end