function plot_SpatialError(u_L1_norm_error_total,v_L1_norm_error_total, u_L2_norm_error_total, v_L2_norm_error_total,L_x,L_y,Lz,N_vec,Test)


switch Test
    case '1D_Adv'
        figure()
        loglog(L_x./N_vec,u_L1_norm_error_total,'o-','MarkerSize',5); hold on; grid on;
        loglog(L_x./N_vec,u_L2_norm_error_total,'o-','MarkerSize',5);
        loglog(L_x./N_vec,(L_x./N_vec.^2),'k--')
        legend({'L1 norm','L2 norm','dh^2'},'Location','northwest')
        ylabel('Error')
        xlabel('Delta grid points (1/N)')
        
    case {'2D_TGV','2D_Vortex'}
        figure()
        subplot(2,1,1)
        loglog(L_x./N_vec,u_L1_norm_error_total,'o-','MarkerSize',5); hold on; grid on;
        loglog(L_x./N_vec,u_L2_norm_error_total,'o-','MarkerSize',5);
        loglog(L_x./N_vec,(L_x./N_vec.^2),'k--')
        legend({'L1 norm','L2 norm','dh^2'},'Location','northwest')
        ylabel('u error')
        xlabel('Delta grid points (1/N)')
        subplot(2,1,2)
        loglog(L_y./N_vec,v_L1_norm_error_total,'o-','MarkerSize',5); hold on; grid on;
        loglog(L_y./N_vec,v_L2_norm_error_total,'o-','MarkerSize',5);
        loglog(L_y./N_vec,(L_y./N_vec.^2),'k--')
        legend({'L1 norm','L2 norm','dh^2'},'Location','northwest')
        ylabel('v error')
        xlabel('Delta grid points (1/N)')

    case '3D_TGV'

end



end
