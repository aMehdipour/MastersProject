% Define the entry dynamics function
function dy = entryDynamics(e, y, sigma0, sigmaf, e_initial, e_final, constants, state, num_steps, CFD)
    % Extract state variables from y
    r = y(1) * constants.LENGTH_SCALE;
    lon = y(2);
    lat = y(3);
    fpa = y(4);
    psi = y(5);
    
    V = sqrt(2 * (1 / r - e));

    % Compute bank angle
    bank = bankAngleProfile(e, sigma0, sigmaf, e_initial, e_final);

    % Compute aerodynamic forces (placeholder for ComputeAero function)
    [L, D, S, ~, ~, ~] = ComputeAero(V, rho, state.aoa, state.ss, flapDeflection, CFD, constants);

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
    V_dot = -D - constants.GRAVITY_NOMINAL * sinfpa / r^2 - constants.OMEGA_EARTH^2 * r * coslat * (sinfpa * coslat - cosfpa * sinlat * cospsi);
    fpa_dot = (L * cosbank) / V - (constants.GRAVITY_NOMINAL * cosfpa) / (V * r) - (2 * constants.OMEGA_EARTH * V * coslat * sinpsi) / (V * r) - constants.OMEGA_EARTH^2 * r * coslat * (cosfpa * coslat + sinfpa * cospsi * sinlat);
    psi_dot = (L * sinbank) / (V * cosfpa) + (V * cosfpa * sinpsi * tan(lat)) / r - (2 * constants.OMEGA_EARTH * V * (tan(fpa) * cospsi * coslat + sinlat)) / (V * cosfpa) + (constants.OMEGA_EARTH^2 * r * sinpsi * sinlat * coslat) / (V * cosfpa);

    % Pack the derivatives into a column vector
    dy = [r_dot; lon_dot; lat_dot; fpa_dot; psi_dot];
end
