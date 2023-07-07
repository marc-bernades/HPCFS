function order = Compute_Spatial_Order_Accuracy(norm_error_total,N_vec)

order = log(norm_error_total(2:end)./norm_error_total(1:end-1))./(log(N_vec(1:end-1)./N_vec(2:end)));

end
