function DataOutput(u,v,w,E,rho,mu,kappa,P,T,ke,e,c_p,c_v,sos,Beta_T,Beta_v,Beta_s,Alpha_p,X,Y,Z,t_vec,ke_total,Invariants,name_file_out,varargin)

% Convert to column
u_data     = u(:);
v_data     = v(:);
w_data     = w(:);
E_data     = E(:);
rho_data   = rho(:);
mu_data    = mu(:);
kappa_data = kappa(:);
P_data     = P(:);
T_data     = T(:);
ke_data    = ke(:);
e_data     = e(:);
c_p_data   = c_p(:);
c_v_data   = c_v(:);
sos_data   = sos(:);
Beta_T_data   = Beta_T(:);
Beta_v_data   = Beta_v(:);
Beta_s_data   = Beta_s(:);
Alpha_p_data  = Alpha_p(:);
X_data     = X(:);
Y_data     = Y(:);
Z_data     = Z(:);

% Save to csv
Headings      = {'u', 'v', 'w', 'E', 'rho', 'mu', 'kappa', 'P', 'T', 'ke', 'e', 'c_p','c_v',...
    'sos', 'Beta_T', 'Beta_v', 'Beta_s', 'Alpha_p','X', 'Y', 'Z'};
Data          = [u_data, v_data, w_data, E_data, rho_data, mu_data, kappa_data,  P_data, T_data, ...
    ke_data, e_data, c_p_data, c_v_data, sos_data, Beta_T_data, Beta_v_data, Beta_s_data, ...
    Alpha_p_data, X_data, Y_data, Z_data];
Data_output   = [Headings; num2cell(Data)];

% Convert cell to a table and use first row as variable names
T = cell2table(Data_output(2:end,:),'VariableNames',Data_output(1,:));
 
% Write the table to a CSV file
writetable(T,[name_file_out + ".csv"])

%% Time variables and Invariants
% Save to csv
Variables_Inv = fieldnames(Invariants);
% for index = 1:length(Variables_Inv)
%     eval([Variables_Inv{index} ' = Invariants. ' Variables_Inv{index} ';']);
% end

% Trim back the pre-allocated values up to the flag
if isempty(varargin)
    b_flag = 0;
else
    b_flag = varargin{1};
    [t_vec,ke_total,Invariants] = TrimTimeVariables(t_vec,ke_total,Invariants,b_flag);
end

% Set data into matrix
Headings      = {'t_vec', 'ke_total'}; Headings = [Headings, Variables_Inv'];
Data          = [t_vec', ke_total'];
for index = 1:length(Variables_Inv)
    Data = [Data, Invariants.(Variables_Inv{index})'];
end
Data_output   = [Headings; num2cell(Data)];

% Convert cell to a table and use first row as variable names
T = cell2table(Data_output(2:end,:),'VariableNames',Data_output(1,:));
 
% Write the table to a CSV file
writetable(T,[name_file_out + "_Time.csv"])


end

