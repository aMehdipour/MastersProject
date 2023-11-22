function dy = predictor(t, y, sigma0, sigmaf, t_initial, t_final, constants)
  % Extract state variables from y
  r = y(1);% * constants.LENGTH_SCALE;
  lon = y(2);
  lat = y(3);
  V = y(4);
  fpa = y(5);
  psi = y(6);

  % V = sqrt(2 * (1 / r - e));

  % Compute bank angle
  bank = 0;%bankAngleProfile(t, sigma0, sigmaf, t_initial, t_final);

  % Compute aerodynamic forces using the dummy model
  [L, D, S] = ComputeAeroDummy(V, constants);
  S = 0;
  L = L / constants.MASS;
  D = D / constants.MASS;

  % Pre-compute trigonometric values
  sinfpa = sin(fpa);
  cosfpa = cos(fpa);
  sinpsi = sin(psi);
  cospsi = cos(psi);
  sinlat = sin(lat);
  coslat = cos(lat);
  cosbank = cos(bank);
  sinbank = sin(bank);

  % EOMs
  r_dot = V * sinfpa;
  lon_dot = (V * cosfpa * sinpsi) / (r * coslat);
  lat_dot = (V * cosfpa * cospsi) / r;
  V_dot = -D - (constants.MU * sinfpa) / r^2 + constants.OMEGA_EARTH^2 * r * coslat * (sinfpa * coslat - cosfpa * sinlat * cospsi);
  fpa_dot = 1 / V * (L * cosbank + (V^2 / r - constants.MU / r^2) * cosfpa + 2 * constants.OMEGA_EARTH * V * coslat * sinpsi);
  psi_dot = 1 / V * ((L * sinbank)/ cosfpa + V^2 / r * cosfpa * sinpsi * tan(lat) - 2 * constants.OMEGA_EARTH * V * (tan(fpa) * cospsi * coslat - sinlat)...
                     + (constants.OMEGA_EARTH^2 * r ) / cosfpa * sinpsi * sinlat * coslat);

  % Pack the derivatives into a column vector
  dy = [r_dot; lon_dot; lat_dot; V_dot; fpa_dot; psi_dot];
end
