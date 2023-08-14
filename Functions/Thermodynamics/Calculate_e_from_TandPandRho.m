function [ e ] = Calculate_e_from_TandPandRho(bSolver, T,rho,P,Fluid, Substance )

if strcmp(bSolver,'Real')

    % PengRobinson
    [a,b,R,dadT,d2adT2,NASA_coefficients] = PengRobinson(T, Substance);

    % Calculate ideal entalpy depending on temperature
    h_ideal = zeros(size(T)); % Preallocate
    for i = 1:length(T(:,1,1))
        for j = 1:length(T(1,:,1))
            for k = 1:length(T(1,1,:))
                if T(i,j,k) >= 200 && T(i,j,k) < 1000
                    h_ideal(i,j,k) = R*T(i,j,k).*(NASA_coefficients(8) + NASA_coefficients(9)*T(i,j,k)/2 + NASA_coefficients(10)*T(i,j,k).^2/3 + NASA_coefficients(11)*T(i,j,k).^3/4 + NASA_coefficients(12)*T(i,j,k).^4/5 + NASA_coefficients(13)./T(i,j,k));
                elseif T(i,j,k) >= 1000 && T(i,j,k) < 6000
                    h_ideal(i,j,k) = R*T(i,j,k).*(NASA_coefficients(1) + NASA_coefficients(2)*T(i,j,k)/2 + NASA_coefficients(3)*T(i,j,k).^2/3 + NASA_coefficients(4)*T(i,j,k).^3/4 + NASA_coefficients(5)*T(i,j,k).^4/5 + NASA_coefficients(6)./T(i,j,k));
                elseif T(i,j,k) < 200
                    h_ideal(i,j,k) = R*200.*(NASA_coefficients(8) + NASA_coefficients(9)*200/2 + NASA_coefficients(10)*200.^2/3 + NASA_coefficients(11)*200.^3/4 + NASA_coefficients(12)*200.^4/5 + NASA_coefficients(13)./200);
                    h_slope        = R*(NASA_coefficients(8) + NASA_coefficients(9)*200 + NASA_coefficients(10)*200.^2 + NASA_coefficients(11)*200.^3 + NASA_coefficients(12)*200.^4);
                    h_ideal(i,j,k) = h_ideal(i,j,k) + h_slope.*(T(i,j,k)-200);
                else
                    break
                end
            end
        end
    end


    % Departure function based on temperature
    % v = 1./rho;
    % K1 = 1/sqrt(8)/b*log((v+(1-sqrt(2))*b)./(v+(1+sqrt(2))*b));
    % dep = (a - T.*dadT).*K1;
    % e = h_ideal - R*T + dep;

    % Departure function based on compressibility factor
    v   = 1./rho;
    Z   = P.*v./(R*T);
    B   = b.*P./(R.*T);
    dep = R.*T.*(Z - 1) + (a - T.*dadT)./(2*sqrt(2)*b).*log((Z + (1 - sqrt(2)).*B)./(Z + (1 + sqrt(2)).*B));
    e   = h_ideal - P.*v + dep;

%     [c_p,c_v,gamma] = Calculate_SpecificHeatCapacities(bSolver, P,rho,T, Fluid,Substance);
%     [sos, Beta_T, Beta_v, Beta_s, Alpha_p] = Calculate_sos(bSolver,rho,T,c_p,Fluid,Substance);
%     Beta_v = 1./T;
%     c_p = 1.0047E+3 + zeros(size(T));
%     e = T./rho.*(rho.*c_p - Beta_v.*P);


else % IDEAL GAS
    e     = P./(rho*(Fluid.gamma - 1));


end
end

