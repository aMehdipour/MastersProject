function sigma = bankAngleProfile(t, sigma0, sigmaf, t_initial, t_final)
    sigma = sigma0 + (sigmaf - sigma0) * (t - t_initial) / (t_final - t_initial);
end
