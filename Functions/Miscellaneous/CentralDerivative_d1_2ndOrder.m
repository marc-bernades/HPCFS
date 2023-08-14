function [dx, dy, dz] = CentralDerivative_d1_2ndOrder(u)
%#codegen


% Allocate memory for derivatives
dx = zeros(length(u(:,1,1)),length(u(1,:,1)),length(u(1,1,:)));
dy = dx;
dz = dx;

%% X direction
% Extremes upwind first order (Boundary points not needed)
dx(:,1,:)    = u(:,2,:) - u(:,1,:);
dx(:,end,:)  = u(:,end,:) - u(:,end - 1,:);
% Sweep internal points
for i = 2:(length(u(1,:,1))-1)
    dx(:,i,:)    = 1/2*(u(:,i+1,:) - u(:,i-1,:));
end

%% Y direction
dy(1,:,:)    = u(2,:,:) - u(1,:,:);
dy(end,:,:)  = u(end,:,:) - u(end - 1,:,:);
for j = 2:(length(u(:,1,1))-1)
    dy(j,:,:)    = 1/2*(u(j+1,:,:) - u(j-1,:,:));
end


%% Z direction
dz(:,:,1)    = u(:,:,2) - u(:,:,1);
dz(:,:,end)  = u(:,:,end) - u(:,:,end - 1);
for k = 2:(length(u(1,1,:))-1)
    dz(:,:,k)    = 1/2*(u(:,:,k+1) - u(:,:,k-1));
end




end

