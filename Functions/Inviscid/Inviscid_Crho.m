function C = Inviscid_Crho(rho,u,v,w,fi,dx,dy,dz)

[dfiu_x,~,~] = CentralDerivative_d1_2ndOrder(fi.*u);
[~,dfiv_y,~] = CentralDerivative_d1_2ndOrder(fi.*v);
[~,~,dfiw_z] = CentralDerivative_d1_2ndOrder(fi.*w);

[drho_x,drho_y,drho_z] = CentralDerivative_d1_2ndOrder(rho);


C = rho.*dfiu_x./(dx)  + rho.*dfiv_y./(dy)  + rho.*dfiw_z./(dz) + ...
        fi.*u.*drho_x./(dx)  + fi.*v.*drho_y./(dy)  + fi.*w.*drho_z./(dz);

% Set outer points to 0
C([1,end],:,:)    = 0;
C(:,[1,end],:)    = 0;
C(:,:,[1,end])    = 0;

end