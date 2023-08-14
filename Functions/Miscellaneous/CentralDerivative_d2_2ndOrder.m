function [dx, dy, dz] = CentralDerivative_d2_2ndOrder(u)


% Allocate memory for derivatives
dx = zeros(length(u(:,1,1)),length(u(1,:,1)),length(u(1,1,:)));
dy = dx;
dz = dx;

%% X direction
% Extremes left to 0 (Boundary points not needed)
% Sweep internal points
for i = 2:(length(u(1,:,1))-1)
    dx(:,i,:)    = (u(:,i+1,:) - 2*u(:,i,:) + u(:,i-1,:));
end

%% Y direction
for j = 2:(length(u(:,1,1))-1)
    dy(j,:,:)    = (u(j+1,:,:) - 2*u(j,:,:) + u(j-1,:,:));
end

%% Z direction
for k = 2:(length(u(1,1,:))-1)
    dz(:,:,k)    = (u(:,:,k+1) - 2*u(:,:,k) + u(:,:,k-1));
end




end

