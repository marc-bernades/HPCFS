function [sos, Beta_T, Beta_v, Beta_s, Alpha_p] = Calculate_sos(bSolver,rho,T,c_p,P,bPressureModel,Fluid,Substance)

if strcmp(bSolver,'Real')

    % PengRobinson
    [a,b,R,dadT,d2adT2,NASA_coefficients] = PengRobinson(T, Substance);

    % Volume
    v = 1./rho;

    % Thermodynamic expressions for sos
    dP_dv_constant_T = -1*R*T./(v - b).^2 + a.*(2*v + 2*b)./(v.^2 + 2*v*b - b^2).^2;
    dP_dT_constant_v = R./(v - b) - dadT./(v.^2 + 2.*v.*b - b.^2);
    Beta_T           = -1./(v.*dP_dv_constant_T); %Isothermal compressability
    Beta_v           = -1*(dP_dT_constant_v./(v.*dP_dv_constant_T)); %Volume expansivity
    Beta_s           = Beta_T - (v.*T.*Beta_v.^2)./c_p; %Isentropic compressability
    sos              = sqrt(1./(rho.*Beta_s));

    Alpha_p          = 1./(v.*dP_dv_constant_T).*dP_dT_constant_v; %Thermal expansivity

    % nasa polynomial for N2
    % an = [3.531005280E+00,-1.236609870E-04,-5.029994370E-07,2.435306120E-09,...
    %      -1.408812350E-12,-1.046976280E+03,2.967474680E+00];
    %
    % v = 1./rho;
    %
    %
    % cp_ideal = R*(an(1) + an(2)*T + an(3)*T.^2 + an(4)*T.^3 + an(5)*T.^4);
    % cv_ideal = cp_ideal - R;
    %
    % dpdT = R./(v-b) - dadT./(v.^2+2*v.*b-b.^2);
    % dpdv = -R.*T./(v-b).^2.*(1-2*a.*((R.*T.*(v+b).*((v.^2+2*v.*b-b.^2)./(v.^2-b.^2)).^2).^(-1)));
    % K1 = 1/sqrt(8)./b.*log((v+(1-sqrt(2)).*b)./(v+(1+sqrt(2)).*b));
    % cp = cv_ideal - T.*(dpdT).^2./dpdv - K1.*T.*d2adT2;
    %
    % kT = -1./v./dpdv;
    % av = -dpdT./v./dpdv;
    % ks = kT - v.*T.*av.^2./cp;
    %
    % sos = sqrt(1./rho./ks);

else
    % Ideal gas

    % Volume
    v = 1./rho;

    % Thermodynamic expressions for sos
    dP_dv_constant_T = -1*Fluid.R_specific*T./(v).^2;
    dP_dT_constant_v = Fluid.R_specific./(v);
    Beta_T           = -1./(v.*dP_dv_constant_T); %Isothermal compressability
    Beta_v           = -1*(dP_dT_constant_v./(v.*dP_dv_constant_T)); %Volume expansivity
    Beta_s           = Beta_T - (v.*T.*Beta_v.^2)./c_p; %Isentropic compressability
    if bPressureModel == 1
        sos          = sqrt(Fluid.gamma*P./rho); % To avoid error when transporting pressure
    else
        sos          = sqrt(1./(rho.*Beta_s)); % idem to sos   = sqrt(Fluid.gamma*p./rho)
    end
    Alpha_p          = 1./(v.*dP_dv_constant_T).*dP_dT_constant_v; %Thermal expansivity


end



end