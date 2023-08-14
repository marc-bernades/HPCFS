function X = stretching_grid(N,x_0,A,L)

% Number of points
num_points = 1:1:N+2;
% eta        = (num_points - 0.5)/N;

% Stretching
% X = x_0 + L*eta + A*(0.5*L - L*eta).*(1 - eta).*eta;

X = zeros(1,length(num_points));

for i = 1:max(num_points)
    eta        = (i - 1 - 0.5)/N;
    X(i) =  x_0 + L*eta + A*(0.5*L - L*eta).*(1 - eta).*eta;

    if i == 1
        eta        = (1 - 0.5)/N;
        X(i) =  x_0 - ( L*eta + A*( 0.5*L - L*eta )*( 1.0 - eta )*eta );


    elseif  i == N + 2
        eta = ( N - 0.5 )/N;
        X(i) = x_0 + 2.0*L - ( L*eta + A*( 0.5*L - L*eta )*( 1.0 - eta )*eta );
    end
end

end