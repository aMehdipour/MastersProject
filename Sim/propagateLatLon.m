function [lat, lon] = propagateLatLon(latOld, lonOld, V, heading, dt, constants)
    % This function propagates forward the vehicle position,
    % stored as latitude and longitude, using the updated
    % velocity and heading. Requires latitude, longitude,
    % and heading to be in radians. Requires velocity to
    % be in m/s.

    % Calculations
    eccSq    = 1 - (constants.RADIUS_POLE^2 / constants.RADIUS_EQ^2); % Earth's eccentricity squared
    N        = constants.RADIUS_EQ / sqrt(1 - eccSq * sin(latOld)^2);
    velNorth = V * cos(heading); % Vehicle North velocity
    velEast  = V * sin(heading); % Vehicle East velocity

    N = constants.RADIUS_EQ / sqrt(1 - eccSq * sin(latOld)^2); % Earth's radius of curvature in the prime vertical
    M = constants.RADIUS_EQ * (1 - eccSq) / (1 - eccSq * sin(latOld)^2)^(3/2); % Earth's radius of curvature in the meridian

    deltaLat = (velNorth / M) * dt;
    deltaLon = (velEast / (N * cos(latOld))) * dt;

    lat = latOld + deltaLat;
    lon = lonOld + deltaLon;