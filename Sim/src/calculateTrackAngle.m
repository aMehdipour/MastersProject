function chi = calculateTrackAngle(V, lat, lon, heading, ss, constants, dt)
  % Calculates the ground track angle, chi. Assumes velocity
  % in m/s and heading in radians.
  eccSq    = 1 - (constants.RADIUS_POLE^2 / constants.RADIUS_EQ^2); % Earth's eccentricity squared
  velNorth = V * cos(heading) * cos(ss) - V * sin(heading) * sin(ss); % Vehicle North velocity
  velEast  = V * sin(heading) * cos(ss) + V * cos(heading) * sin(ss); % Vehicle East velocity

  N = constants.RADIUS_EQ / sqrt(1 - eccSq * sin(lat)^2); % Earth's radius of curvature in the prime vertical

  % Account for Earth's rotation in the North velocity
  velNorthRelECEF = velNorth + constants.OMEGA_EARTH * N * cos(lat);

  chi = atan2(velEast, velNorthRelECEF);
