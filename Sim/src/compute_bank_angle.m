function sigma = compute_bank_angle(sigma_0, sigma_f, e, e_0, e_f)
  sigma = sigma_0 + (sigma_f - sigma_0) * (e - e_0) / (e_f - e_0);
