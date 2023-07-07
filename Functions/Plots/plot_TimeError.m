function plot_TimeError(u_L1_norm_error_total,v_L1_norm_error_total, u_L2_norm_error_total, v_L2_norm_error_total,dt_vec,RK_order,Test)

switch Test
    case '1D_Adv'
        figure()
        loglog(dt_vec,u_L1_norm_error_total,'o-','MarkerSize',5); hold on; grid on;
        loglog(dt_vec,u_L2_norm_error_total,'o-','MarkerSize',5);
        loglog(dt_vec,(dt_vec.^RK_order),'k--')
        legend({'L1 norm','L2 norm',["dt^" + num2str(RK_order)]},'Location','northwest')
        ylabel('rho error')
        xlabel('Time step (s)')
        
    case '2D_TGV'
        figure()
        subplot(2,1,1)
        loglog(dt_vec,u_L1_norm_error_total,'o-','MarkerSize',5); hold on; grid on;
        loglog(dt_vec,u_L2_norm_error_total,'o-','MarkerSize',5);
        loglog(dt_vec,(dt_vec.^RK_order),'k--')
        legend({'L1 norm','L2 norm',["dt^" + num2str(RK_order)]},'Location','northwest')
        ylabel('u error')
        xlabel('Time step (s)')
        subplot(2,1,2)
        loglog(dt_vec,v_L1_norm_error_total,'o-','MarkerSize',5); hold on; grid on;
        loglog(dt_vec,v_L2_norm_error_total,'o-','MarkerSize',5);
        loglog(dt_vec,(dt_vec.^RK_order),'k--')
        legend({'L1 norm','L2 norm',["dt^" + num2str(RK_order)]},'Location','northwest')
        ylabel('v error')
        xlabel('Time step (s)')
        
    case '3D_TGV'

end




end
