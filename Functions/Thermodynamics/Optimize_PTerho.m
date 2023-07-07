function [P,T] = Optimize_PTerho(P_current,T_current,e_target,rho_target,Fluid,Substance)

% Initial guess x0 as current state
x = [P_current,T_current];
x0  = x;

%% Non-linear fsolve Method
% Define optimisation function
    function F = root2d(x)
        F    = zeros(2,1);
        F(1) = Calculate_Rho_from_TandP('Real',x(2),x(1),Fluid,Substance) - rho_target;
        F(2) = Calculate_e_from_TandPandRho('Real',x(2),F(1)+rho_target,x(1), Fluid,Substance) - e_target;
%         F(2) = Calculate_e_from_TandPandRho_Test('Real',x(2),F(1)+rho_target,x(1), Fluid,Substance,c_p,1/x(2)) - e_target;
    end

% Call solver
fun     = @root2d;
options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-8,'MaxIterations',100); %,'PlotFcn',@optimplotfirstorderopt);
% options = optimoptions('fsolve','algorithm','levenberg-marquardt','Display','off','FunctionTolerance',1e-8,'MaxIterations',100); %,'PlotFcn',@optimplotfirstorderopt);
x       = fsolve(fun,x0,options);

%% Newton Raphson Method
% fun = @(x) [Calculate_Rho_from_TandP(x(2),x(1),Substance) - rho_target, Calculate_e_from_TandPandRho(x(2),rho_target,x(1), Substance) - e_target];
% % fprintf('\ninitial guess: P = %g[Pa], T = %g[K]\n',x0) % display initial guess 
% options = optimset('TolX',1e-8); % set TolX
% [x, resnorm, f, exitflag, output, jacob] = newtonraphson(fun, x0, options);

%% Assign results
P = x(1);
T = x(2);

end
