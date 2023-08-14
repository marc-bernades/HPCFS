function [u_L1_norm_error, v_L1_norm_error] = Calculate_Error_L1(u,v,u_exact,v_exact,hx,hy,hz,Test)

% Remove outer points
% u_exact = u_exact(2:end-1,2:end-1,2:end-1);
% v_exact = v_exact(2:end-1,2:end-1,2:end-1);
% u       = u(2:end-1,2:end-1,2:end-1);
% v       = v(2:end-1,2:end-1,2:end-1);

% u_L1_norm_error = sum(sum(sum(abs(u_exact - u))))*(hx*hy*hz)./sum(sum(sum(abs(u_exact))));
% v_L1_norm_error = sum(sum(sum(abs(v_exact - v))))*(hx*hy*hz)./sum(sum(sum(abs(v_exact))));

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

u_L1_norm_error = sum(sum(sum(abs(u_exact - u))))*(h_total);
v_L1_norm_error = sum(sum(sum(abs(v_exact - v))))*(h_total);

end

