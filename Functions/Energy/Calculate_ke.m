function ke_total = Calculate_ke(u, v, w, X, Y, Z)
    
    % Geometric dimensions
    [dx,~,~] = CentralDerivative_d1_2ndOrder(X);
    [~,dy,~] = CentralDerivative_d1_2ndOrder(Y);
    [~,~,dz] = CentralDerivative_d1_2ndOrder(Z);

    % Volume
    volume = dx.*dy.*dz;

    % kinetic energy
    ke = 0.5*(u.*u + v.*v + w.*w);
    ke_volume = ke.*volume;

    % Sum across the domain
    ke_volume_sum = sum(sum(sum(ke_volume(2:end-1,2:end-1,2:end-1))));
    volume_sum    = sum(sum(sum(volume(2:end-1,2:end-1,2:end-1))));

    % Total energy
    ke_total = ke_volume_sum/volume_sum;


end

