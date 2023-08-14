function [ke, e, P, T, c_p, c_v, gamma, sos, Beta_T, Beta_v, Beta_s, Alpha_p,E] = Update_ThermodynamicState(bSolver,rho,u,v,w,E,P,T,Fluid,Substance,bPressureModel)

ke  = 0.5*(u.^2 + v.^2 + w.^2);
% E   = ke + e; % TEC, non-PEP: Uncomment this, comment other row and swap input e for E
e   = E - ke;

if bPressureModel == 1
    % Pressure obtained already > Hence decoupled from T.
    % We perform Aiteken with T
    T     = Calculate_T_fromPandRho(bSolver,P,rho,Fluid,Substance);
    e     = Calculate_e_from_TandPandRho(bSolver, T,rho,P,Fluid, Substance );
    E     = e + ke;

else
    % Update T and P with internal energy and EoS
    if strcmp(bSolver,'Real')

        % Implement optimisation algorithm > Iterate point by point
        for i = 1:length(u(:,1,1))
            for j = 1:length(u(1,:,1))
                for k = 1:length(u(1,1,:))
                    % First guess as previous step of P and T
                    [P(i,j,k),T(i,j,k)] = Optimize_PTerho(P(i,j,k),T(i,j,k),e(i,j,k),rho(i,j,k),Fluid,Substance);
                end
            end
        end


    else
        % For ideal gas
        P = e.*rho*(Fluid.gamma - 1);
        T = Calculate_T_fromPandRho(bSolver,P,rho,Fluid,Substance);
        %  T = e./c_v;

    end

end

[c_p,c_v,gamma]                          = Calculate_SpecificHeatCapacities(bSolver,P,rho,T, Fluid,Substance);
[sos, Beta_T, Beta_v, Beta_s, Alpha_p]   = Calculate_sos(bSolver,rho,T,c_p,P,bPressureModel,Fluid,Substance);

end

