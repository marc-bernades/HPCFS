function order = Compute_Temp_Order_Accuracy(norm_error_total,dt_vec)

order = log(norm_error_total(2:end)./norm_error_total(1:end-1))./(log((dt_vec(2:end))./dt_vec(1:end-1)));

end
