function [c_p,c_v,gamma] = Calculate_SpecificHeatCapacities(bSolver, P,rho,T, Fluid,Substance)

if strcmp(bSolver,'Real')
    % PengRobinson
    [a,b,R,dadT,d2adT2,NASA_coefficients] = PengRobinson(T, Substance);

    %% Cp real gas
    % Calculate ideal Cp depending on temperature
    c_p_ideal = zeros(size(T)); % Preallocate
    for i = 1:length(T(:,1,1))
        for j = 1:length(T(1,:,1))
            for k = 1:length(T(1,1,:))
                if T(i,j,k) >= 200 && T(i,j,k) < 1000
                    c_p_ideal(i,j,k) = R.*(NASA_coefficients(8) + NASA_coefficients(9)*T(i,j,k) + NASA_coefficients(10)*T(i,j,k).^2 + NASA_coefficients(11)*T(i,j,k).^3 + NASA_coefficients(12)*T(i,j,k).^4);
                elseif T(i,j,k) >= 1000 && T(i,j,k) < 6000
                    c_p_ideal(i,j,k) = R.*(NASA_coefficients(1) + NASA_coefficients(2)*T(i,j,k) + NASA_coefficients(3)*T(i,j,k).^2 + NASA_coefficients(4)*T(i,j,k).^3 + NASA_coefficients(5)*T(i,j,k).^4);
                elseif T(i,j,k) < 200
                    % Assume constant temperature below 200K
                    c_p_ideal(i,j,k) = R.*(NASA_coefficients(8) + NASA_coefficients(9)*200 + NASA_coefficients(10)*200.^2 + NASA_coefficients(11)*200.^3 + NASA_coefficients(12)*200.^4);
                else
                    break
                end
            end
        end
    end

    % Departure function Cp
    v         = 1./rho;
    Z         = P.*v./(R*T);
    A         = a.*P./(R.*T).^2;
    B         = b.*P./(R.*T);
    M         = (Z.^2 + 2*B.*Z - B.^2)./(Z - B);
    N         = dadT.*(B./(b*R));
    dep_c_p   = (R.*(M - N).^2)./(M.^2 - 2*A.*(Z + B)) - (T.*d2adT2./(2*sqrt(2)*b)).*log((Z + (1 - sqrt(2))*B)./(Z + (1 + sqrt(2))*B)) - R;

    % Cp real gas
    c_p       = c_p_ideal + dep_c_p;


    % an = [3.531005280E+00,-1.236609870E-04,-5.029994370E-07,2.435306120E-09,...
    %      -1.408812350E-12,-1.046976280E+03,2.967474680E+00];
    %
    % v = 1./rho;
    %
    % T = getTfromPandRho(p,rho);
    % [a,b,R,dadT,d2adT2] = getThermo(T);
    % cp_ideal = R*(an(1) + an(2)*T + an(3)*T.^2 + an(4)*T.^3 + an(5)*T.^4);
    % cv_ideal = cp_ideal - R;
    %
    % dpdT = R./(v-b) - dadT./(v.^2+2*v.*b-b.^2);
    % dpdv = -R.*T./(v-b).^2.*(1-2*a.*((R.*T.*(v+b).*((v.^2+2*v.*b-b.^2)./(v.^2-b.^2)).^2).^(-1)));
    % K1 = 1/sqrt(8)./b.*log((v+(1-sqrt(2)).*b)./(v+(1+sqrt(2)).*b));
    % cp = cv_ideal - T.*(dpdT).^2./dpdv - K1.*T.*d2adT2;

    %% Cv real gas

    % Cv ideal gas
    c_v_ideal = c_p_ideal - R;

    % Departure function Cv
    dep_c_v = -1*(T.*d2adT2)./(2+sqrt(2).*b).*log((Z + (1 - sqrt(2))*B)./(Z + (1 + sqrt(2))*B));

    % Cv real gas
    c_v  = c_v_ideal + dep_c_v;

else

    % Ideal gas
    c_v   = zeros(size(P)) + Fluid.R_specific/(Fluid.gamma - 1);
    c_p   = zeros(size(P)) + c_v*Fluid.gamma;
end

%% Gamma
gamma = c_p./c_v;

end