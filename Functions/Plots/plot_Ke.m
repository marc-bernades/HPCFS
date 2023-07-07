function plot_Ke(t_vec,ke_total,varargin)

if isempty(varargin)
    figure()
    plot(t_vec, ke_total/ke_total(1))
    xlabel('time [s]')
    ylabel('ke/ke_0')
    title("Kinetic Energy Evolution")

else
    figure()
    plot(t_vec, ke_total); hold on
    plot(varargin{1},varargin{2})
    xlabel('time [s]')
    ylabel('ke/ke_0')
    title("Kinetic Energy Evolution")
    legend('KGP','Shima')

end

end