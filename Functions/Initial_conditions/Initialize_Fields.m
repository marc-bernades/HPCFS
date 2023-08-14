function [u,v,w,P,T] = Initialize_Fields(bSolver,X,Y,Z,L_x,L_y,L_z,Test,Fluid,Substance)

switch Test

%     case '1D_Adv'
%         %% 1D Adv test
%         u = 0*X + Fluid.U_0;
%         v = 0*u;
%         w = 0*u;
%         rho = (Fluid.rho_min + Fluid.rho_max)/2  + (Fluid.rho_max - Fluid.rho_min)/2*sin(2*pi*X);
%         P = Fluid.P_0+0*X;
%         T = P./(rho*Fluid.R_specific);

    case '1D_HighPressure'
        %% 1D High Pressure
        u = 0*X;
        v = 0*u;
        w = 0*u;
        P = Fluid.P_0 + 0*X;
%         T = Fluid.T_0*X + 100;
        T = Fluid.T_0*X + Fluid.T_0;
        % Iterate T if rho is input
        % T = Calculate_T_fromPandRho( p,rho,Substance)

    case '1D_Adv_Real'
        %% 1D Adv non-linear ref Ma
        [rho, P, T, e, u,v,w ke, E ] = Reference_Solution_1D_Adv(bSolver,X,Fluid,Substance);

    case '1D_Adv_Shima'
        %% 1D Adv non-linear ref Ma
        [rho, P, T, e, u,v,w ke, E ] = Reference_Solution_1D_Adv_Shima(bSolver,X,Fluid,Substance);

    case '1D_Adv'
        %% 1D Adv non-linear ref Ma
        [rho, P, T, e, u,v,w ke, E ] = Reference_Solution_1D_Adv_Shima(bSolver,X,Fluid,Substance);


    case 'Lid-Driven_Cavity'
        %% Lid-driven cavity
        u = 0*X;
        v = 0*u;
        w = 0*u;
        P = Fluid.P_0+0*X;
        T = P./(Fluid.rho_0*Fluid.R_specific);


    case '2D_Vortex'

        %% 2D Vortex centered at 0,0 with radius 1. Assume Mach v = March inf
        r = sqrt((X - Fluid.Vortex_x0).^2 + (Y - Fluid.Vortex_y0).^2);
        r_norm = r/Fluid.r_v;
        u = Fluid.U_0.*(1 - Fluid.Ma_vortex./Fluid.Ma*(((Y - Fluid.Vortex_y0)./Fluid.r_v).*exp((1 - r_norm.^2)/2)));
        v = Fluid.U_0.*(Fluid.Ma_vortex./Fluid.Ma*(((X - Fluid.Vortex_x0)./Fluid.r_v).*exp((1 - r_norm.^2)/2)));
        w = u*0;
        P = Fluid.P_0.*((1 - ((Fluid.gamma - 1)/2*Fluid.Ma_vortex^2).*exp((1 - r_norm.^2))).^(Fluid.gamma/(Fluid.gamma - 1)));
        rho = Fluid.rho_0*((1 - ((Fluid.gamma - 1)/2*Fluid.Ma_vortex^2).*exp((1 - r_norm.^2))).^(1/(Fluid.gamma - 1)));
        T = P./(rho*Fluid.R_specific);
        %[~, ~, ~, T, ~] = Update_ThermodynamicState(rho,u,v,w,E,Fluid.gamma,Fluid.c_v);

    case '2D_MixingLayer'

        %% 2D Mixing Layer reference
        % Velocity slip
        delta_u     = 0.20;     % Velocity amplitude
        delta_T     = 0.75/2;   % T amplitude to reach from 0.75xTc to 1.5Tc
        a           = 20;
        u = Fluid.U_0.*(1 + delta_u*tanh(a*Y)); % u = zeros(size(X)) - 0.5*tanh(4*pi*Y) or u = u + A_0.*(rand(size(X)).*2 -1);
        v = u*0;
        w = u*0;

        % Perturbation on velocity
        [~,dy,~]        = CentralDerivative_d1_2ndOrder(Y);
        [u_pert,v_pert] = Calculate_2D_MixingLayer_Perturbation(u,dy,X,Y,L_y);
        % Set pertrubation y limits
        y_limit = 0.1; % To tanh non-linearity
        [~,Index_top] = min(abs(Y(:,2,2) - y_limit));
        [~,Index_bot] = min(abs(Y(:,2,2) + y_limit));
        u(Index_bot:Index_top,:,:) = u(Index_bot:Index_top,:,:) + 1*u_pert(Index_bot:Index_top,:,:);
        v(Index_bot:Index_top,:,:) = v(Index_bot:Index_top,:,:) + 1*v_pert(Index_bot:Index_top,:,:);
        
        % Thermodynamics
        P   = zeros(size(X)) + Fluid.P_0;
        [~, T_c,  ~, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Substance);
        T = T_c*(0.75 + delta_T - delta_T*tanh(a*Y));
        % Rho input
%         rho = rho_c*(1 + delta_rho*tanh(a*Y));
        % Set temperature to get Rspecific
%         T_getThermo = 300;
%         [~,~,R,~,~,~] = PengRobinson(T_getThermo, Substance);
        % Estimate first guess
%         vol = 1./rho;
%         T = P.*vol./R; % Estimation ideal gas
%         T = Calculate_T_fromPandRho( bSolver, P,rho,T,Fluid,Substance); % T   = zeros(size(X)) +  200 + 200*tanh(4*pi*Y);

    case '2D_MixingLayer_Periodic'

        %% 2D Mixing Layer reference
        % Velocity slip
        delta_T     = 0.75/2;   % T amplitude to reach from 0.75xTc to 1.5Tc
        a           = 20;
        u = zeros(size(Y));
        delta = pi/15;
        u(Y<=pi) = tanh((Y(Y<=pi) - 0.5*pi)/delta);
        u(Y>pi)  = tanh((3/2*pi - Y(Y>pi))/delta);
        Epsilon = 0.05;
        v = Epsilon*sin(X);
        w = u*0;

        % Thermodynamics
        P   = zeros(size(X)) + Fluid.P_0;
        T   = zeros(size(X)) + Fluid.T_0;
%         [~, T_c,  ~, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Substance);
%         T = T_c*(0.75 + delta_T - delta_T*tanh(a*Y));

    case '2D_MixingLayer_Periodic_Real'

        %% 2D Mixing Layer reference
        % Velocity slip
        u = zeros(size(Y));
        delta = pi/15;
        u(Y<=pi) = tanh((Y(Y<=pi) - 0.5*pi)/delta);
        u(Y>pi)  = tanh((3/2*pi - Y(Y>pi))/delta);
        Epsilon = 0.05;
        v = Epsilon*sin(X);
        w = u*0;

        % Thermodynamics
        P   = zeros(size(X)) + Fluid.P_0;
        T   = zeros(size(X));
        [~, T_c,  ~, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Substance);
        % Change from 0.75 to 1.5Tc
        T(Y<=pi) = 1.125 + 0.75*0.5*tanh((Y(Y<=pi) - 0.5*pi)/delta);
        T(Y>pi)  = 1.125 + 0.75*0.5*tanh((3/2*pi - Y(Y>pi))/delta);
        T = T*T_c;
    
    case '2D_TGV'

        %% 2D TGV reference
        u = Fluid.U_0*cos(X).*sin(Y);
        v = -Fluid.U_0*sin(X).*cos(Y);
        w = u*0;
        P = Fluid.P_0 + Fluid.rho_0.*(Fluid.U_0*Fluid.U_0/4).*(cos(2*X) + cos(2*Y));
        T = P./(Fluid.rho_0*Fluid.R_specific);

    case '2D_ChannelFlow'

        % 2D Channel flow
        n_random = 2*rand(size(X)) - 1;
        y_dist = min(Y, 2*Fluid.delta - Y);
        u = 2*Fluid.u_0*y_dist/Fluid.delta + Fluid.alpha*Fluid.u_0*n_random;
        v = u*0;
        w = u*0;

        % Thermodynamics
        P   = zeros(size(X)) + Fluid.P_0;
        T = Fluid.T_bw + Y./(2*Fluid.delta)*(Fluid.T_tw - Fluid.T_bw);

    case '3D_TGV'
        %% 3D TGV JCP
        u = sin(X).*cos(Y).*cos(Z);
        v = -cos(X).*sin(Y).*cos(Z);
        w = u*0;
        P = Fluid.P_0 + (1/16).*((cos(2*X) + cos(2*Y)).*(cos(2*Z) + 2) - 2);
        T = P./(Fluid.rho_0.*Fluid.R_specific);


        %     case '3D_TGV'
        %         %% 3D TGV reference
        %         u = U_0*(2/sqrt(3))*sin(2*pi/3)*sin(X).*cos(Y).*cos(Z);
        %         v = U_0*(2/sqrt(3))*sin(-2*pi/3)*cos(X).*sin(Y).*cos(Z);
        %         w = U_0*(2/sqrt(3))*sin(0)*cos(X).*cos(Y).*sin(Z);
        %         P = P_0 + rho_0.*(U_0*U_0/16).*(cos(2*X) + cos(2*Y)).*(cos(2*Z) + 2);
        %         T = P./(rho_0.*R_specific);

    case '3D_TGV_Real'
        %% 3D TGV JCP
        u = sin(X).*cos(Y).*cos(Z);
        v = -cos(X).*sin(Y).*cos(Z);
        w = u*0;
        P = Fluid.P_0 + (1/16).*((cos(2*X) + cos(2*Y)).*(cos(2*Z) + 2) - 2);
        T = Fluid.T_0 + 0*u;
%         rho = Fluid.rho_0 + 0*P;
%         T   = Calculate_T_fromPandRho( bSolver,P,rho,Fluid,Substance);
    case '3D_TGV_Viscous'
        %% AMR 2019 - Brachet JFM 1983
        theta = 0;
        u = Fluid.U_0*2/sqrt(3)*sin(theta + 2/3*pi).*sin(X).*cos(Y).*cos(Z);
        v = Fluid.U_0*2/sqrt(3)*sin(theta - 2/3*pi).*cos(X).*sin(Y).*cos(Z);
        w = Fluid.U_0*2/sqrt(3)*sin(theta).*cos(X).*cos(Y).*cos(Z);
        P = Fluid.P_0 + (1/16).*((cos(2*X) + cos(2*Y)).*(cos(2*Z) + 2) - 2);
        T = P./(Fluid.rho_0.*Fluid.R_specific);

end


end