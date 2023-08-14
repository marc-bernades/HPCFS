function C = Inviscid_CL(rho,u,v,w,fi,dx,dy,dz)

[du_x,~,~] = CentralDerivative_d1_2ndOrder(u);
[~,dv_y,~] = CentralDerivative_d1_2ndOrder(v);
[~,~,dw_z] = CentralDerivative_d1_2ndOrder(w);

[drho_x,drho_y,drho_z] = CentralDerivative_d1_2ndOrder(rho);

[dfi_x,dfi_y,dfi_z] = CentralDerivative_d1_2ndOrder(fi);

C  = rho.*fi.*du_x./(dx)  + rho.*fi.*dv_y./(dy)  + rho.*fi.*dw_z./(dz) + ...
        rho.*u.*dfi_x./(dx)  + rho.*v.*dfi_y./(dy)  + rho.*w.*dfi_z./(dz) + ...
        fi.*u.*drho_x./(dx)  + fi.*v.*drho_y./(dy)  + fi.*w.*drho_z./(dz);

% Set outer points to 0
C([1,end],:,:)    = 0;
C(:,[1,end],:)    = 0;
C(:,:,[1,end])    = 0;

end