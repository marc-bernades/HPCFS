function Conv = Calculate_Convective_Split(C_D, C_u, C_fi, C_rho, C_L, coeff)

alpha   = coeff(1);
beta    = coeff(2);
gamma   = coeff(3);
delta   = coeff(4);
epsilon = coeff(5);


Conv    = alpha*C_D + beta*C_fi + gamma*C_u + delta*C_rho + epsilon*C_L;


end

