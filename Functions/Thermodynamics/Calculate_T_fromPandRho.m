function T = Calculate_T_fromPandRho( bSolver, p,rho,Fluid,Substance)

if strcmp(bSolver,'Real')

    v = 1./rho;

    % Set temperature to get Rspecific
    T_getThermo = 300;
    [~,~,R,~,~,~] = PengRobinson(T_getThermo, Substance);

    % Estimate first as ideal gas
    T = p.*v./R;

    %% Aitken optimizer
    % Tolerance and max iterations
    tolerance        = 1E-8;
    n_iterations_max = 1000;

    for i = 1:length(T(:,1,1))
        for j = 1:length(T(1,:,1))
            for k = 1:length(T(1,1,:))
                % Aitken iteration
                for n = 1:n_iterations_max
                    x0 = T(i,j,k);
                    [a,b,R,~,~,~] = PengRobinson(x0, Substance);
                    % Peng Robinson equation
                    x1 = (v(i,j,k) - b)/R.*(p(i,j,k) + a/(v(i,j,k)^2 + 2*b*v(i,j,k) - b^2));
                    [a,b,R,~,~,~] = PengRobinson(x1, Substance);
                    x2 = (v(i,j,k) - b)/R.*(p(i,j,k) + a/(v(i,j,k)^2 + 2*b*v(i,j,k) - b^2));
                    T(i,j,k) = x2 - ((x2 - x1)^2)/(x2 - 2*x1 + x0);
                    if abs((T(i,j,k) - x2)/T(i,j,k)) < tolerance
                        break
                    end

                end

            end
        end
    end


%     T = p.*v./Fluid.R_specific;

else
    % Ideal gas
    v = 1./rho;
    T = p.*v./Fluid.R_specific;

end
end
