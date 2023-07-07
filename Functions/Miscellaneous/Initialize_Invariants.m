function[t_vec, ke_total, Invariants] = Initialize_Invariants(b_flag,margin,t_end,t_0,dt)

t_vec                  = zeros(1,round((t_end-t_0)/dt)*margin) + b_flag;
ke_total               = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rho_norm    = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rhou_bar    = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rhoE_norm   = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rhoe_norm   = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rhos_norm   = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rhoui2_norm = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rhoE2_norm  = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rhoe2_norm  = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.rhos2_norm  = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.T_norm      = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.e_norm      = zeros(1,round((t_end-t_0)/dt)*margin);
Invariants.P_norm      = zeros(1,round((t_end-t_0)/dt)*margin);


end