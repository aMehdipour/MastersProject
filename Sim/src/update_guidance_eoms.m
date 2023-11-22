function guidance_state = update_guidance_eoms(state, guidance_state, CFD, rho, flapDeflection, num_steps, constants)
  % Unpack what we need from the state vector
  r               = guidance_state.r * constants.LENGTH_SCALE;
  lon             = guidance_state.longitude;
  lat             = guidance_state.latitude;
  fpa             = guidance_state.flightPathAngle;
  psi             = guidance_state.heading;
  bank            = guidance_state.bank;

  de = 1 / num_steps;
  
  % Calculate aero forces based on given state
  V = sqrt(2 * (1 / r - guidance_state.e));
  [L, D, S, ~, ~, ~] = ComputeAero(V, rho, state.aoa, state.ss, flapDeflection, CFD, constants);

  % Pre-compute some values
  sinfpa  = sin(fpa);
  cosfpa  = cos(fpa);
  sinpsi  = sin(psi);
  cospsi  = cos(psi);
  sinlat  = sin(lat);
  coslat  = cos(lat);
  cosbank = cos(bank);
  sinbank = sin(bank);

  % Equations of motion
  r_dot     = V * sinfpa;
  lon_dot   = (V * cosfpa * sinpsi) / (r * coslat);
  lat_dot   = (V * cosfpa * cospsi) / r;
  V_dot     = -D - constants.GRAVITY_NOMINAL * sinfpa / r^2 - constants.OMEGA_EARTH^2 * r * coslat * (sinfpa * coslat - cosfpa * sinlat * cospsi);
  fpa_dot   = (L * cosbank) / V - (constants.GRAVITY_NOMINAL * cosfpa) / (V * r) - (2 * constants.OMEGA_EARTH * V * coslat * sinpsi) / (V * r) - constants.OMEGA_EARTH^2 * r * coslat * (cosfpa * coslat + sinfpa * cospsi * sinlat);
  psi_dot   = (L * sinbank) / (V * cosfpa) + (V * cosfpa * sinpsi * tan(lat)) / r - (2 * constants.OMEGA_EARTH * V * (tan(fpa) * cospsi * coslat + sinlat)) / (V * cosfpa) + (constants.OMEGA_EARTH^2 * r * sinpsi * sinlat * coslat) / (V * cosfpa);
  s_dot = -V * cosfpa / r;

  guidance_state.r   = r + r_dot * de;
  guidance_state.lon = lon + lon_dot * de;
  guidance_state.lat = lat + lat_dot * de;
  guidance_state.V   = V + V_dot * de;
  guidance_state.fpa = fpa + fpa_dot * de;
  guidance_state.psi = psi + psi_dot * de;
  guidance_state.s   = s + s_dot * de;
  guidance_state.e   = 1 / guidance_state.r - guidance_state.V^2 / 2;

