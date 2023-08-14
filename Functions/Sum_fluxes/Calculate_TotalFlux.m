function [rho_tot, rhou_tot, rhov_tot, rhow_tot, rhoE_tot] = Calculate_TotalFlux(rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv, ...
                                                                      dP_rhou, dP_rhov, dP_rhow, dP_rhoE,...
                                                                      rhou_visc, rhov_visc, rhow_visc, rhoE_visc,bEnergySplit, bEnergySplitEnthalpy, bFlux, ...
                                                                      f_rhou, f_rhov, f_rhow, f_rhoE, f_rhouvw)

rho_tot  = -rho_conv;
if bFlux == 1
    % Pressure within convective in flux form
    rhou_tot = -rhou_conv + rhou_visc;
    rhov_tot = -rhov_conv + rhov_visc;
    rhow_tot = -rhow_conv  + rhow_visc;

else
    rhou_tot = -rhou_conv - dP_rhou + rhou_visc;
    rhov_tot = -rhov_conv - dP_rhov + rhov_visc;
    rhow_tot = -rhow_conv - dP_rhow + rhow_visc;
end
if bEnergySplit == 1 ||  bEnergySplitEnthalpy == 1 || bFlux == 1
    % rhoE includes pressure gradient
    rhoE_tot = -rhoE_conv + rhoE_visc;
else
    rhoE_tot = -rhoE_conv - dP_rhoE + rhoE_visc;
end

% Add source terms update only inner points
rhou_tot(2:end-1,2:end-1,2:end-1) = rhou_tot(2:end-1,2:end-1,2:end-1) + f_rhou(2:end-1,2:end-1,2:end-1);
rhov_tot(2:end-1,2:end-1,2:end-1) = rhov_tot(2:end-1,2:end-1,2:end-1) + f_rhov(2:end-1,2:end-1,2:end-1);
rhow_tot(2:end-1,2:end-1,2:end-1) = rhow_tot(2:end-1,2:end-1,2:end-1) + f_rhow(2:end-1,2:end-1,2:end-1);
rhoE_tot(2:end-1,2:end-1,2:end-1) = rhoE_tot(2:end-1,2:end-1,2:end-1) + f_rhoE(2:end-1,2:end-1,2:end-1) + f_rhouvw(2:end-1,2:end-1,2:end-1); 



end

