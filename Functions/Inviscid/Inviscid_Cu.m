function C = Inviscid_Cu(rho,u,v,w,fi,dx,dy,dz)

[du_x,~,~] = CentralDerivative_d1_2ndOrder(u);
[~,dv_y,~] = CentralDerivative_d1_2ndOrder(v);
[~,~,dw_z] = CentralDerivative_d1_2ndOrder(w);

[drhofi_x,drhofi_y,drhofi_z] = CentralDerivative_d1_2ndOrder(rho.*fi);


C = u.*drhofi_x./(dx)  + v.*drhofi_y./(dy)  + w.*drhofi_z./(dz) + ...
         rho.*fi.*(du_x./(dx) + dv_y./(dy) + dw_z./(dz));
 

% Set outer points to 0
C([1,end],:,:)    = 0;
C(:,[1,end],:)    = 0;
C(:,:,[1,end])    = 0;

end