function [guidance_state, error] = predictor_step(state, sigma_0, e_0, e_f, g, rho, CFD, flapDeflection, constants)

  guidance_state = state;
  num_steps = ceil((e_f - e_0) / dt);

  e = e_0;

  for step = 1:num_steps
    guidance_state.bank = compute_bank_angle(sigma_0, e, e_0, e_f, constants.SIGMA_F);
    guidance_state = update_guidance_eoms(guidance_state, e, CFD, rho, flapDeflection, g, dt, constants);
  end

  error = 0.5 * (guidance_state.s - constants.RANGE_FINAL)^2;  
