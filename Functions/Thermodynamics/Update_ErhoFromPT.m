function [E,rho] = Update_ErhoFromPT(u,v,w,P,T,gamma,c_v)
    e   = T.*c_v; %(1/rho*(gamma-1))*P
    rho = P./(e.*(gamma - 1)); % (1/(R_specific*T))*P
    ke  = 0.5*(u.^2 + v.^2 + w.^2);
    E   = e + ke;
end