function [F, J] = funwrapper(fun, x)
% if nargout<2 use finite differences to estimate J
try
    [F, J] = fun(x);
catch
    F = fun(x);
    J = jacobian(fun, x); % evaluate center diff if no Jacobian
end
F = F(:); % needs to be a column vector
end
