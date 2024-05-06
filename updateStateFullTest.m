function [state, derivatives] = updateStateFull(state, commandedBankAngle, constants, parameters)
    % Extracting state variables
    velocity = state.velocity;
    flightPathAngle = state.flightPathAngle;
    bankAngle = state.bankAngle;
    latitude = state.latitude;
    longitude = state.longitude;
    heading = state.heading;
    normGeocentricDistance = state.normGeocentricDistance;
    rangeToGo = state.rangeToGo;
    dt = normalizeState(parameters.dt, 'time', true, constants);
    velocityUnnormalized = state.velocityUnnormalized;
    altitudeUnnormalized = state.altitudeUnnormalized;
    earthRotationRate = normalizeState(constants.EARTH_ROTATION_RATE, 'time', true, constants);

    % Precalculate trigonometric functions
    sin_flightPathAngle = sin(flightPathAngle);
    cos_flightPathAngle = cos(flightPathAngle);
    tan_flightPathAngle = tan(flightPathAngle);
    sin_latitude = sin(latitude);
    cos_latitude = cos(latitude);
    tan_latitude = tan(latitude);
    cos_longitude = cos(longitude);
    sin_heading = sin(heading);
    cos_heading = cos(heading);
    sin_bankAngle = sin(bankAngle);
    cos_bankAngle = cos(bankAngle);

    % Compute non-dimensionalized lift and drag forces
    [rho, ~, ~] = atmosphereModel(altitudeUnnormalized);
    [~, trimCD, trimCL] = calculateTrimAero(parameters, constants, velocityUnnormalized, altitudeUnnormalized);
    lift = 0.5 * rho * trimCL * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;
    drag = 0.5 * rho * trimCD * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;

    % Calculate heading error relative to the terminal point
    [~, terminalAzimuth] = distance('gc', latitude, longitude, deg2rad(constants.TARGET_LAT), deg2rad(constants.TARGET_LON), 'radians');
    state.headingError = terminalAzimuth - heading;

    % Equations of motion using precomputed trigonometric values
    derivatives.velocityDot = -drag - sin_flightPathAngle / normGeocentricDistance^2 + earthRotationRate^2 * normGeocentricDistance * cos_latitude * (sin_flightPathAngle * cos_latitude - cos_flightPathAngle * sin_latitude * cos_heading);
    derivatives.flightPathAngleDot = (lift * cos_bankAngle + (velocity^2 - 1 / normGeocentricDistance) * cos_flightPathAngle / normGeocentricDistance + 2 * earthRotationRate * velocity * cos_latitude * sin_heading + earthRotationRate^2 * normGeocentricDistance * cos_latitude * (cos_flightPathAngle * cos_latitude + sin_flightPathAngle * cos_heading * sin_latitude)) / velocity;
    derivatives.normGeocentricDistanceDot = velocity * sin_flightPathAngle;
    derivatives.rangeDot = -velocity * cos_flightPathAngle / normGeocentricDistance;
    derivatives.lonDot = (velocity * cos_flightPathAngle * sin_heading) / (normGeocentricDistance * cos_longitude);
    derivatives.latDot = (velocity * cos_flightPathAngle * cos_heading) / normGeocentricDistance;
    derivatives.headingDot = 1 / velocity * ((lift * sin_bankAngle) / cos_flightPathAngle + (velocity^2 * cos_flightPathAngle * sin_heading * tan_latitude) / normGeocentricDistance - 2 * earthRotationRate * velocity * (tan_flightPathAngle * cos_heading * cos_latitude - sin_latitude) + (earthRotationRate^2 * normGeocentricDistance * sin_heading * sin_latitude * cos_latitude) / cos_flightPathAngle);
    derivatives.bankAngleDot = clampBankAngleRate(bankAngle, commandedBankAngle, constants, parameters.dt);

    % Update state variables
    state.velocity = velocity + derivatives.velocityDot * dt;
    state.flightPathAngle = flightPathAngle + derivatives.flightPathAngleDot * dt;
    state.normGeocentricDistance = normGeocentricDistance + derivatives.normGeocentricDistanceDot * dt;
    state.rangeToGo = rangeToGo + derivatives.rangeDot * dt;
    state.bankAngle = bankAngle + derivatives.bankAngleDot;
    state.latitude = latitude + derivatives.latDot * dt;
    state.longitude = longitude + derivatives.lonDot * dt;
    state.heading = heading + derivatives.headingDot * dt;
    state.lift = lift;
    state.drag = drag;
    state.loadFactor = sqrt(lift^2 + drag^2);
    state.heatingRate = constants.K_Q * sqrt(rho) * velocity^3.15 * 1e-4;
    state.posEcef = lla2ecef([state.latitude, state.longitude, state.altitudeUnnormalized], 'WGS84');
    state.altitudeUnnormalized = normalizeState(state.normGeocentricDistance,'length',false,constants);
    state.velocityUnnormalized = normalizeState(state.velocity,'velocity',false,constants);
    state.rangeToGoUnnormalized = normalizeState(state.rangeToGo, 'length', false, constants);
    state.rotNedFromEcef = calcRotNedFromEcef(state.latitude, state.longitude);
    state.posNed = state.rotNedFromEcef * state.posEcef';
    state.crossrange = calculateCrossrange(state, constants);
    state.deadbandWidth = calculateDeadbandWidth(state.velocity, constants);
end
