function [u_exact,v_exact] = Calculate_Exact_2D_TGV(X,Y,Z,time,mu,rho,U_0)

nu = mu./rho;

u_exact =  U_0*cos(X).*sin(Y).*exp(-2*nu.*time);
v_exact = -U_0*sin(X).*cos(Y).*exp(-2*nu.*time);

% u_exact = U_0*(2/sqrt(3))*sin(2*pi/3)*sin(X).*cos(Y).*cos(Z).*exp(-2*nu.*time);
% v_exact = U_0*(2/sqrt(3))*sin(-2*pi/3)*cos(X).*sin(Y).*cos(Z).*exp(-2*nu.*time);


end

