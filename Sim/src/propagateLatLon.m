function [lat, lon] = propagateLatLon(latOld, lonOld, V, heading, dt, constants)
    % This function propagates forward the vehicle position,
    % stored as latitude and longitude, using the updated
    % velocity and heading.
    % Requires latitude, longitude, and heading to be in radians.
    % Requires velocity to be in m/s.

    % FIXME: This function is not accurate for large time steps. Also add in altitude as an input

    % Calculations
    eccSq    = 1 - (constants.RADIUS_POLE^2 / constants.RADIUS_EQ^2); % Earth's eccentricity squared
    velNorth = V * cos(heading); % Vehicle North velocity
    velEast  = V * sin(heading); % Vehicle East velocity

    N = constants.RADIUS_EQ / sqrt(1 - eccSq * sin(latOld)^2); % Earth's radius of curvature in the prime vertical
    M = constants.RADIUS_EQ * (1 - eccSq) / (1 - eccSq * sin(latOld)^2)^(3/2); % Earth's radius of curvature in the meridian

    deltaLat = (velNorth / (M + altitude)) * dt;
    deltaLon = (velEast / ((N + altitude) * cos(latOld))) * dt;

    lat = latOld + deltaLat;
    lon = lonOld + deltaLon;

    % Wrap the longitude to the [-180, 180) degree interval.
    lon = mod(lon + pi, 2 * pi) - pi;

    % Limit the latitude to the [-90, 90] degree interval.
    lat = max(min(lat, pi / 2), -pi / 2);