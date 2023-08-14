function Y = TimeIntegration_RKGeneral(Y_0,dt,k,i_RK,RK_order,a,b,c)


if i_RK < RK_order
    % Term of aij·F(u1)
    a_F = 0;
    for ii = 1:i_RK
        a_F = a_F + a(i_RK+1,ii)*k{ii};
    end
    % Update of ui_RK+1
    Y = Y_0 + dt*a_F;

else
    %Update final step of RK with b coefficients bjj·kjj
    b_k = 0;
    for jj = 1:length(b)
        b_k = b_k + b(jj)*k{jj};
    end
    Y = Y_0 + dt*b_k;



end

end