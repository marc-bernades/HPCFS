function plot_TurbulentSpectra(k, e_k)


figure
loglog(k,e_k,'linewidth',2)
hold on
loglog(k,k.^(-5/3),'--')
legend('KGP','Intertial range')
xlabel('${{\omega}_{k}}$','interpreter','latex')
ylabel('${{k_e}({\omega}_{k})}$','interpreter','latex')

grid on


end