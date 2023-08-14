function [rho_conv, rhou_conv, rhov_conv, rhow_conv, rhoE_conv] = Calculate_Convective(C_D, C_u, C_fi, C_rho, C_L, scheme)

switch scheme
    case 0 %Div
        rho_conv  = C_D{1};
        rhou_conv = C_D{2};
        rhov_conv = C_D{3};
        rhow_conv = C_D{4};
        rhoE_conv = C_D{5};

    case 1 %KGP
        rho_conv  = 1/2*(C_D{1} + C_u{1});
        rhou_conv = 1/4*(C_D{2} + C_fi{2} + C_u{2} + C_rho{2});
        rhov_conv = 1/4*(C_D{3} + C_fi{3} + C_u{3} + C_rho{3});
        rhow_conv = 1/4*(C_D{4} + C_fi{4} + C_u{4} + C_rho{4});
        rhoE_conv = 1/4*(C_D{5} + C_fi{5} + C_u{5} + C_rho{5});

    case 2 % Feireisen F form (KCP 2019)

        % Weights Feiereisen F form
        zi    = 1;
        delta = 0;

        % Coefficient constraints
        alpha   = 1/2 - delta;
        beta    = zi/2;
        gamma   = delta;
        epsilon = (1 - zi) / 2 - delta;

        % Convective terms
        rho_conv  = zi*C_D{1} + (1 - zi)*C_u{1};
        rhou_conv = alpha*C_D{2} + beta*C_fi{2} + gamma*C_u{2} + delta*C_rho{2} + epsilon*C_L{2};
        rhov_conv = alpha*C_D{3} + beta*C_fi{3} + gamma*C_u{3} + delta*C_rho{3} + epsilon*C_L{3};
        rhow_conv = alpha*C_D{4} + beta*C_fi{4} + gamma*C_u{4} + delta*C_rho{4} + epsilon*C_L{4};
        rhoE_conv = alpha*C_D{5} + beta*C_fi{5} + gamma*C_u{5} + delta*C_rho{5} + epsilon*C_L{5};
    
    case 3 % C

        % Weights Feiereisen F form
        zi    = 0;
        delta = 1/2;

        % Coefficient constraints
        alpha   = 1/2 - delta;
        beta    = zi/2;
        gamma   = delta;
        epsilon = (1 - zi) / 2 - delta;

        % Convective terms
        rho_conv  = zi*C_D{1} + (1 - zi)*C_u{1};
        rhou_conv = alpha*C_D{2} + beta*C_fi{2} + gamma*C_u{2} + delta*C_rho{2} + epsilon*C_L{2};
        rhov_conv = alpha*C_D{3} + beta*C_fi{3} + gamma*C_u{3} + delta*C_rho{3} + epsilon*C_L{3};
        rhow_conv = alpha*C_D{4} + beta*C_fi{4} + gamma*C_u{4} + delta*C_rho{4} + epsilon*C_L{4};
        rhoE_conv = alpha*C_D{5} + beta*C_fi{5} + gamma*C_u{5} + delta*C_rho{5} + epsilon*C_L{5};

    otherwise
        rho_conv  = C_D{1};
        rhou_conv = C_D{2};
        rhov_conv = C_D{3};
        rhow_conv = C_D{4};
        rhoE_conv = C_D{5};
end

% Set outer points to 0
rho_conv([1,end],:,:)    = 0;
rho_conv(:,[1,end],:)    = 0;
rho_conv(:,:,[1,end])    = 0;
rhou_conv([1,end],:,:)   = 0;
rhou_conv(:,[1,end],:)   = 0;
rhou_conv(:,:,[1,end])   = 0;
rhov_conv([1,end],:,:)   = 0;
rhov_conv(:,[1,end],:)   = 0;
rhov_conv(:,:,[1,end])   = 0;
rhow_conv([1,end],:,:)   = 0;
rhow_conv(:,[1,end],:)   = 0;
rhow_conv(:,:,[1,end])   = 0;
rhoE_conv([1,end],:,:)   = 0;
rhoE_conv(:,[1,end],:)   = 0;
rhoE_conv(:,:,[1,end])   = 0;





end

