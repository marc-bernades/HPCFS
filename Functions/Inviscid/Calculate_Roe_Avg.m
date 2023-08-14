function U = Calculate_Roe_Avg(u,dx,dimension)

% Compute variable at +1/2
if dimension == 1
    u_2  = circshift(u,[0,1,0]);
    dx_2 = circshift(dx,[0,1,0]);
    U    = (u.*dx_2 + u_2.*dx)./(dx + dx_2);

elseif dimension == 2
    u_2  = circshift(u,[1,0,0]);
    dx_2 = circshift(dx,[1,0,0]);
    U    = (u.*dx_2 + u_2.*dx)./(dx + dx_2);

elseif dimension == 3
    u_2  = circshift(u,[1,0,0]);
    dx_2 = circshift(dx,[0,0,1]);
    U    = (u.*dx_2 + u_2.*dx)./(dx + dx_2);

end



end