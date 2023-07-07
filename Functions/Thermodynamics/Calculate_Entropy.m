function s = Calculate_Entropy(T,rho,P,Substance)

% Call PengRobinson
[~,b,R,dadT,~,NASA_coefficients] = PengRobinson(T, Substance);
    

%% S real gas
    % Calculate ideal S depending on temperature
    s_ideal = zeros(size(T)); % Preallocate
    for i = 1:length(T(:,1,1))
        for j = 1:length(T(1,:,1))
            for k = 1:length(T(1,1,:))
                if T(i,j,k) >= 200 && T(i,j,k) < 1000
                    s_ideal(i,j,k) = R.*(NASA_coefficients(8)*log(T(i,j,k)) + NASA_coefficients(9)*T(i,j,k) + NASA_coefficients(10)*1/2*T(i,j,k).^2 + NASA_coefficients(11)*1/3*T(i,j,k).^3 + NASA_coefficients(12)*1/4*T(i,j,k).^4 + NASA_coefficients(14));
                elseif T(i,j,k) >= 1000 && T(i,j,k) < 6000
                    s_ideal(i,j,k) = R.*(NASA_coefficients(1)*log(T(i,j,k)) + NASA_coefficients(2)*T(i,j,k) + NASA_coefficients(3)*1/2*T(i,j,k).^2 + NASA_coefficients(4)*1/3*T(i,j,k).^3 + NASA_coefficients(5)*1/4*T(i,j,k).^4 + NASA_coefficients(7));
                elseif T(i,j,k) < 200
                    % Assume constant temperature below 200K
                    s_ideal(i,j,k) = R.*(NASA_coefficients(8)*log(T(i,j,k)) + NASA_coefficients(9)*T(i,j,k) + NASA_coefficients(10)*1/2*T(i,j,k).^2 + NASA_coefficients(11)*1/3*T(i,j,k).^3 + NASA_coefficients(12)*1/4*T(i,j,k).^4 + NASA_coefficients(14));
                else
                    break
                end
            end
        end
    end

    % Departure function s ideal
    v          = 1./rho;
    Z          = P.*v./(R*T);
    B          = b.*P./(R.*T);
    dep_s      = R*log(Z - B) + P.*dadT./(2*sqrt(2)*R*T.*B).*log((Z + 2.414*B)./(Z - 2.414*B));

    % S real gas
    s       = s_ideal + dep_s;



end