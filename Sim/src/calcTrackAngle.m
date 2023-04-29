function chi = calcTrackAngle(V, lat, lon, heading, constants)
    %   Calculates the track angle of the vehicle based on the current
    %   velocity vector, latitude, longitude, and heading.
    velNorth = V * cos(heading);
    velEast  = V * sin(heading);
