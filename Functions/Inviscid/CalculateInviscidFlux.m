function [rho_inv_flux,rhou_inv_flux,rhov_inv_flux,rhow_inv_flux,rhoE_inv_flux] = CalculateInviscidFlux(rho,u,v,w,E,sos,P,dx,dy,dz,bFlux,bSolver, Fluid, Substance)

% Flux scheme
scheme  = bFlux.type;

% Initialize
rho_inv_flux  = zeros(size(u));
rhou_inv_flux = zeros(size(u));
rhov_inv_flux = zeros(size(u));
rhow_inv_flux = zeros(size(u));
rhoE_inv_flux = zeros(size(u));


% Only inner points
for i = 2:length(u(:,1,1))-1
    for j = 2:length(u(1,:,1))-1
        for k = 2:length(u(1,1,:))-1
            %% X - Direction (j component)
            % x = j+1/2
            index_L = j; index_R = j + 1;
            rho_L = rho(i,index_L,k); rho_R = rho(i,index_R,k);
            u_L   = u(i,index_L,k);   u_R   = u(i,index_R,k);
            v_L   = v(i,index_L,k);   v_R   = v(i,index_R,k);
            w_L   = w(i,index_L,k);   w_R   = w(i,index_R,k);
            E_L   = E(i,index_L,k);   E_R   = E(i,index_R,k);
            P_L   = P(i,index_L,k);   P_R   = P(i,index_R,k);
            a_L   = sos(i,index_L,k); a_R   = sos(i,index_R,k);
            
            % rho
            var_type = 0;
            rho_F_p  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhou
            var_type = 1;
            rhou_F_p  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhov
            var_type = 2;
            rhov_F_p  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhow
            var_type = 3;
            rhow_F_p  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhoE
            var_type = 4;
            rhoE_F_p  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);

            % x = j-1/2
            index_L = j - 1; index_R = j;
            rho_L = rho(i,index_L,k); rho_R = rho(i,index_R,k);
            u_L   = u(i,index_L,k);   u_R   = u(i,index_R,k);
            v_L   = v(i,index_L,k);   v_R   = v(i,index_R,k);
            w_L   = w(i,index_L,k);   w_R   = w(i,index_R,k);
            E_L   = E(i,index_L,k);   E_R   = E(i,index_R,k);
            P_L   = P(i,index_L,k);   P_R   = P(i,index_R,k);
            a_L   = sos(i,index_L,k); a_R   = sos(i,index_R,k);
            
            % rho
            var_type = 0;
            rho_F_m  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhou
            var_type = 1;
            rhou_F_m  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhov
            var_type = 2;
            rhov_F_m  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhow
            var_type = 3;
            rhow_F_m  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhoE
            var_type = 4;
            rhoE_F_m  = Riemann_solver(rho_L,rho_R,u_L,u_R,v_L,v_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);

            % Compute flux x-direction
            rho_inv_flux(i,j,k)  = (rho_F_p - rho_F_m)./dx(i,j,k);
            rhou_inv_flux(i,j,k) = (rhou_F_p - rhou_F_m)./dx(i,j,k);
            rhov_inv_flux(i,j,k) = (rhov_F_p - rhov_F_m)./dx(i,j,k);
            rhow_inv_flux(i,j,k) = (rhow_F_p - rhow_F_m)./dx(i,j,k);
            rhoE_inv_flux(i,j,k) = (rhoE_F_p - rhoE_F_m)./dx(i,j,k);

            %% Y - Direction
            % y = i+1/2
            index_L = i; index_R = i + 1;
            rho_L = rho(index_L,j,k); rho_R = rho(index_R,j,k);
            u_L   = u(index_L,j,k);   u_R   = u(index_R,j,k);
            v_L   = v(index_L,j,k);   v_R   = v(index_R,j,k);
            w_L   = w(index_L,j,k);   w_R   = w(index_R,j,k);
            E_L   = E(index_L,j,k);   E_R   = E(index_R,j,k);
            P_L   = P(index_L,j,k);   P_R   = P(index_R,j,k);
            a_L   = sos(index_L,j,k); a_R   = sos(index_R,j,k);
           
            % rho
            var_type = 0;
            rho_F_p  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhou
            var_type = 2;
            rhou_F_p  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhov
            var_type = 1;
            rhov_F_p  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhow
            var_type = 3;
            rhow_F_p  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhoE
            var_type = 4;
            rhoE_F_p  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);

            % y = i-1/2
            index_L = i-1; index_R = i;
            rho_L = rho(index_L,j,k); rho_R = rho(index_R,j,k);
            u_L   = u(index_L,j,k);   u_R   = u(index_R,j,k);
            v_L   = v(index_L,j,k);   v_R   = v(index_R,j,k);
            w_L   = w(index_L,j,k);   w_R   = w(index_R,j,k);
            E_L   = E(index_L,j,k);   E_R   = E(index_R,j,k);
            P_L   = P(index_L,j,k);   P_R   = P(index_R,j,k);
            a_L   = sos(index_L,j,k); a_R   = sos(index_R,j,k);
            
            % rho
            var_type = 0;
            rho_F_m  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhou
            var_type = 2;
            rhou_F_m  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhov
            var_type = 1;
            rhov_F_m  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhow
            var_type = 3;
            rhow_F_m  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhoE
            var_type = 4;
            rhoE_F_m  = Riemann_solver(rho_L,rho_R,v_L,v_R,u_L,u_R,w_L,w_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);

            % Compute flux x-direction
            rho_inv_flux(i,j,k)  = rho_inv_flux(i,j,k)  + (rho_F_p - rho_F_m)./dy(i,j,k);
            rhou_inv_flux(i,j,k) = rhou_inv_flux(i,j,k) + (rhou_F_p - rhou_F_m)./dy(i,j,k);
            rhov_inv_flux(i,j,k) = rhov_inv_flux(i,j,k) + (rhov_F_p - rhov_F_m)./dy(i,j,k);
            rhow_inv_flux(i,j,k) = rhow_inv_flux(i,j,k) + (rhow_F_p - rhow_F_m)./dy(i,j,k);
            rhoE_inv_flux(i,j,k) = rhoE_inv_flux(i,j,k) + (rhoE_F_p - rhoE_F_m)./dy(i,j,k);

            %% Z - Direction
            % z = k+1/2
            index_L = k; index_R = k + 1;
            rho_L = rho(i,j,index_L); rho_R = rho(i,j,index_R);
            u_L   = u(i,j,index_L);   u_R   = u(i,j,index_R);
            v_L   = v(i,j,index_L);   v_R   = v(i,j,index_R);
            w_L   = w(i,j,index_L);   w_R   = w(i,j,index_R);
            E_L   = E(i,j,index_L);   E_R   = E(i,j,index_R);
            P_L   = P(i,j,index_L);   P_R   = P(i,j,index_R);
            a_L   = sos(i,j,index_L); a_R   = sos(i,j,index_R);
            % rho
            var_type = 0;
            rho_F_p  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhou
            var_type = 3;
            rhou_F_p  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhov
            var_type = 2;
            rhov_F_p  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhow
            var_type = 1;
            rhow_F_p  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhoE
            var_type = 4;
            rhoE_F_p  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);

            % z = k-1/2
            index_L = k - 1; index_R = k;
            rho_L = rho(i,j,index_L); rho_R = rho(i,j,index_R);
            u_L   = u(i,j,index_L);   u_R   = u(i,j,index_R);
            v_L   = v(i,j,index_L);   v_R   = v(i,j,index_R);
            w_L   = w(i,j,index_L);   w_R   = w(i,j,index_R);
            E_L   = E(i,j,index_L);   E_R   = E(i,j,index_R);
            P_L   = P(i,j,index_L);   P_R   = P(i,j,index_R);
            a_L   = sos(i,j,index_L); a_R   = sos(i,j,index_R);
            % rho
            var_type = 0;
            rho_F_m  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhou
            var_type = 3;
            rhou_F_m  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhov
            var_type = 2;
            rhov_F_m  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhow
            var_type = 1;
            rhow_F_m  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);
            % rhoE
            var_type = 4;
            rhoE_F_m  = Riemann_solver(rho_L,rho_R,w_L,w_R,v_L,v_R,u_L,u_R,E_L,E_R,P_L,P_R,a_L,a_R,var_type,scheme, bSolver, Fluid, Substance);

            % Compute flux x-direction
            rho_inv_flux(i,j,k)  = rho_inv_flux(i,j,k)  + (rho_F_p - rho_F_m)./dz(i,j,k);
            rhou_inv_flux(i,j,k) = rhou_inv_flux(i,j,k) + (rhou_F_p - rhou_F_m)./dz(i,j,k);
            rhov_inv_flux(i,j,k) = rhov_inv_flux(i,j,k) + (rhov_F_p - rhov_F_m)./dz(i,j,k);
            rhow_inv_flux(i,j,k) = rhow_inv_flux(i,j,k) + (rhow_F_p - rhow_F_m)./dz(i,j,k);
            rhoE_inv_flux(i,j,k) = rhoE_inv_flux(i,j,k) + (rhoE_F_p - rhoE_F_m)./dz(i,j,k);


        end
    end
end


end