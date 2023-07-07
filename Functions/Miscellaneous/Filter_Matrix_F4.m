function    [A,B] = Filter_Matrix_F4(U)

%% Second order filter (Reference Visbal 2002 JFM)
alpha_f = 0.495;
a0      =  5/8  + 3/4*alpha_f;
a1      =  1/2  + alpha_f;
a2      = -1/8  + 1/4*alpha_f;

% Inatialize matricies each direction
Dimension = size(U);
for d = 1:length(Dimension)
    A{d} = zeros([Dimension(d) Dimension(d)]);
    B{d} = zeros([Dimension(d) Dimension(d)]);
end

% Set filter matrices for each direction
for d = 1:length(Dimension)
    if Dimension(d) <= 5 + 2
        A{d} = eye(Dimension(d));
        B{d} = eye(Dimension(d));
    else
        % Set boundary points to prevent singular matrix
        A{d}(1,1) = 1;  A{d}(end,end) = 1;
        B{d}(1,1) = 1;  B{d}(end,end) = 1;

        % Sweep internal points
        for j = 3:Dimension(d)-2
            A{d}(j,j) = 1;
            A{d}(j,j+1) = alpha_f;
            A{d}(j,j-1) = alpha_f;

            B{d}(j,j)   = a0;
            B{d}(j,j+1) = a1/2;
            B{d}(j,j-1) = a1/2;
        end

        % Periodic boundary
        A{d}(2,2) = 1;         A{d}(2,3) = alpha_f;         A{d}(2,end-1) = alpha_f;
        A{d}(end-1,end-1) = 1; A{d}(end-1,end-2) = alpha_f; A{d}(end-1,2) = alpha_f;

        % B matrix d direction
        v = zeros(1,Dimension(d)); v(1) = a2/2; v(2) = a1/2; v(3) = a0; v(4) = a1/2; v(5) = a2/2;
        D = gallery('circul',v);
        B{d}(4:end-3,:) = D(2:end-5,:);
        B{d}(2,end-2) = a2/2; B{d}(2,end-1) = a1/2; B{d}(2,2) = a0; B{d}(2,3) = a1/2; B{d}(2,4) = a2/2;
        B{d}(3,end-1) = a2/2; B{d}(3,2) = a1/2; B{d}(3,3) = a0; B{d}(3,4) = a1/2; B{d}(3,5) = a2/2;
        B{d}(end-2,end-4) = a2/2; B{d}(end-2,end-3) = a1/2; B{d}(end-2,end-2) = a0; B{d}(end-2,end-1) = a1/2; B{d}(end-2,2) = a2/2;
        B{d}(end-1,end-3) = a2/2; B{d}(end-1,end-2) = a1/2; B{d}(end-1,end-1) = a0; B{d}(end-1,2) = a1/2; B{d}(end-1,3) = a2/2;


    end
end

end