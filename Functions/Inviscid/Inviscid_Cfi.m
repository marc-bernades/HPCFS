function C = Inviscid_Cfi(rho,u,v,w,fi,dx,dy,dz)


[drhou_x,~,~] = CentralDerivative_d1_2ndOrder(rho.*u);
[~,drhov_y,~] = CentralDerivative_d1_2ndOrder(rho.*v);
[~,~,drhow_z] = CentralDerivative_d1_2ndOrder(rho.*w);

[dfi_x,dfi_y,dfi_z] = CentralDerivative_d1_2ndOrder(fi);


C = fi.*drhou_x./(dx)  + fi.*drhov_y./(dy)  + fi.*drhow_z./(dz) + ...
        rho.*u.*dfi_x./(dx)  + rho.*v.*dfi_y./(dy)  + rho.*w.*dfi_z./(dz);

% Set outer points to 0
C([1,end],:,:)    = 0;
C(:,[1,end],:)    = 0;
C(:,:,[1,end])    = 0;

end