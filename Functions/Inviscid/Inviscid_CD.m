function C = Inviscid_CD(rho,u,v,w,fi,dx,dy,dz)

[drhoufi_x,~,~] = CentralDerivative_d1_2ndOrder(rho.*u.*fi);
[~,drhovfi_y,~] = CentralDerivative_d1_2ndOrder(rho.*v.*fi);
[~,~,drhowfi_z] = CentralDerivative_d1_2ndOrder(rho.*w.*fi);

C = drhoufi_x./dx + drhovfi_y./dy + drhowfi_z./dz;

% Set outer points to 0
C([1,end],:,:)    = 0;
C(:,[1,end],:)    = 0;
C(:,:,[1,end])    = 0;

end