function [u_L2_norm_error, v_L2_norm_error] = Calculate_Error_L2(u,v,u_exact,v_exact,hx,hy,hz,Test)

% Remove outer points
% u_exact = u_exact(2:end-1,2:end-1,2:end-1);
% v_exact = v_exact(2:end-1,2:end-1,2:end-1);
% u       = u(2:end-1,2:end-1,2:end-1);
% v       = v(2:end-1,2:end-1,2:end-1);

u_L2_norm_error = sum(sum(sum((u_exact - u).^2)));
v_L2_norm_error = sum(sum(sum((v_exact - v).^2)));

%u_L2_norm_den = sum(sum(sum((u_exact).^2)));
%v_L2_norm_den = sum(sum(sum((v_exact).^2)));

%u_L2_norm_error = sqrt((hx*hy*hz)*u_L2_norm_error/u_L2_norm_den);
%v_L2_norm_error = sqrt((hx*hy*hz)*v_L2_norm_error/v_L2_norm_den);
switch Test
    case {'1D_Adv','2D_Vortex'}
        % Grid space
        h_total = hx;
    case '2D_TGV'
        % Area integral
        h_total = hx*hy;
    case '3D_TGV'
        % Volume integral
        h_total = hx*hy*hz;
    otherwise
        h_total = hx*hy*hz;
end

u_L2_norm_error = sqrt((h_total)*u_L2_norm_error);
v_L2_norm_error = sqrt((h_total)*v_L2_norm_error);
end

