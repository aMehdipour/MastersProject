function posECEF = geodeticToECEF(lat, lon, alt, constants)
    a = constants.EARTH_RADIUS_EQ;
    b = constants.EARTH_RADIUS_POLE;
    f = constants.FLATTENING_FACTOR;
    h = alt;

    posECEF = [(a / sqrt(1 - (2*f - f^2) * sin(lat)^2)) * cos(lat) * cos(lon) + h * cos(lat) * cos(lon); ...
               (a / sqrt(1 - (2*f - f^2) * sin(lat)^2)) * cos(lat) * sin(lon) + h * cos(lat) * sin(lon); ...
               (a * (1 - f)^2 / sqrt(1 - (2*f - f^2) * sin(lat)^2)) * sin(lat) + h * sin(lat)];


