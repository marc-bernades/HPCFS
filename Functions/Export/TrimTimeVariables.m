function [t_vec,ke_total,Invariants] = TrimTimeVariables(t_vec,ke_total,Invariants,b_flag)

%% Trim back the pre-allocated values up to the flag
[~,d] = min(abs(t_vec(2:end) - b_flag));
t_vec    = t_vec(1:d); t_vec = t_vec - t_vec(1);
ke_total = ke_total(1:d);
Variables_Inv = fieldnames(Invariants);
for ii = 1:length(Variables_Inv)
    Invariants.(Variables_Inv{ii}) = Invariants.(Variables_Inv{ii})(1:d);
end

end