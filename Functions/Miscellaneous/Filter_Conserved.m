function    [U_filt] = Filter_Conserved(U,A,B,dx,dy,dz,Filter_Type)

% Initialize filtered quantites
U_filt  = U;

A_i = A{1};
A_j = A{2};
A_k = A{3};
B_i = B{1};
B_j = B{2};
B_k = B{3};

switch Filter_Type
    case 'Implicit'
        %% Update filtered quantities
        % For Y direction (i-th)
        for j = 1:length(U(1,:,1))
            for k = 1:length(U(1,1,:))
                filt_temp           = A_i\B_i*U(:,j,k);
                U_filt(2:end-1,j,k) = filt_temp(2:end-1);
            end
        end
        % For X direction (j-th)
        for i = 1:length(U(:,1,1))
            for k = 1:length(U(1,1,:))
                filt_temp           = A_j\B_j*U_filt(i,:,k)';
                U_filt(i,2:end-1,k) = filt_temp(2:end-1)';
            end
        end

        % For Z direction (k-th)
        for i = 1:length(U(:,1,1))
            for j = 1:length(U(1,:,1))
                filt_temp           = A_k\B_k*reshape(U_filt(i,j,:),[length(U(i,j,:)) 1]);
                U_filt(i,j,2:end-1) = filt_temp(2:end-1);
            end
        end


    case 'Gaussian'

        %% Gaussian filter
        [d2u_x,d2u_y,d2u_z]  = CentralDerivative_d2_2ndOrder(U);
        Lap_U                = d2u_x./dx.^2 + d2u_y./dy.^2 + d2u_z./dz.^2;
        grid_size            = (dx.*dy.*dz).^(1/3);
        % Only inner points
        U_filt(2:end-1,2:end-1,2:end-1) = U(2:end-1,2:end-1,2:end-1) + (grid_size(2:end-1,2:end-1,2:end-1).^2)/24.*Lap_U(2:end-1,2:end-1,2:end-1);


end


end
