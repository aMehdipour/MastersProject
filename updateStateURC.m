% Update state variables function
function [state, derivatives] = updateState(state, commandedBankAngle, constants, parameters)
    velocity               = state.velocity;
    flightPathAngle        = state.flightPathAngle;
    angleOfAttack          = state.angleOfAttack;
    sideslipAngle          = state.sideslipAngle;
    latitude               = state.latitude;
    longitude              = state.longitude;
    heading                = state.heading;
    normGeocentricDistance = state.normGeocentricDistance;
    rangeToGo              = state.rangeToGo;
    dt                     = normalizeState(parameters.dt, 'time', true, constants);
    velocityUnnormalized   = state.velocityUnnormalized;
    altitudeUnnormalized   = state.altitudeUnnormalized;
    % crossrange             = state.crossrange;

    % Compute non-dimensionalized lift and drag forces
    [rho, ~, ~] = atmosphereModel(altitudeUnnormalized);
    [trimCL, trimCD, trimCS] = calculateAeroCoefficients(parameters, angleOfAttack, sideslipAngle);
    lift = 0.5 * rho * trimCL * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;
    drag = 0.5 * rho * trimCD * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;

     % Calculate heading error relative to the terminal point
    [~, terminalAzimuth] = distance('gc', latitude, longitude, deg2rad(constants.TARGET_LAT), deg2rad(constants.TARGET_LON), 'radians');
    state.headingError = terminalAzimuth - heading;

    % Equations of motion
    derivatives.velocityDot               = -drag - sin(flightPathAngle) / normGeocentricDistance^2;
    derivatives.flightPathAngleDot        = (lift * cos(commandedBankAngle) + (velocity^2 - 1 / normGeocentricDistance) * cos(flightPathAngle) / normGeocentricDistance) / velocity;
    derivatives.normGeocentricDistanceDot = velocity * sin(flightPathAngle);
    derivatives.rangeDot                  = -velocity * cos(flightPathAngle) / normGeocentricDistance;
    derivatives.lonDot                    = (velocity * cos(flightPathAngle) * sin(heading)) / (normGeocentricDistance * cos(longitude));
    derivatives.latDot                    = (velocity * cos(flightPathAngle) * cos(heading)) / normGeocentricDistance;
    derivatives.headingDot                = 1 / velocity * ((lift * sin(bankAngle)) / cos(flightPathAngle) + (velocity^2) / normGeocentricDistance * cos(flightPathAngle) * sin(heading) * tan(latitude));
    derivatives.bankAngleDot              = clampBankAngleRate(bankAngle, commandedBankAngle, constants, parameters.dt);
    % derivatives.crossrangeDot             = velocity * cos(flightPathAngle) * sin(state.headingError) / normGeocentricDistance;


    % Update state variables
    state.velocity               = velocity + derivatives.velocityDot * dt;
    state.flightPathAngle        = flightPathAngle + derivatives.flightPathAngleDot * dt;
    state.normGeocentricDistance = normGeocentricDistance + derivatives.normGeocentricDistanceDot * dt;
    state.rangeToGo              = rangeToGo + derivatives.rangeDot * dt;
    state.bankAngle              = bankAngle + derivatives.bankAngleDot;
    state.latitude               = latitude + derivatives.latDot * dt;
    state.longitude              = longitude + derivatives.lonDot * dt;
    state.heading                = heading + derivatives.headingDot * dt;
    state.lift                   = lift;
    state.drag                   = drag;
    state.loadFactor             = sqrt(lift^2 + drag^2);
    state.heatingRate            = constants.K_Q * sqrt(rho) * velocity^3.15 * 1e-4;
    state.posEcef                = lla2ecef([state.latitude, state.longitude, state.altitudeUnnormalized], 'WGS84');
    state.altitudeUnnormalized   = normalizeState(state.normGeocentricDistance,'length',false,constants) - constants.EARTH_RADIUS_EQ;
    state.velocityUnnormalized   = normalizeState(state.velocity,'velocity',false,constants);
    state.rangeToGoUnnormalized  = normalizeState(state.rangeToGo, 'length', false, constants);
    state.rotNedFromEcef         = calcRotNedFromEcef(state.latitude, state.longitude);
    state.posNed                 = state.rotNedFromEcef * state.posEcef';
    % [~, state.crossrange]        = calculateDownrangeCrossrange(state,constants);
    % state.crossrange             = normalizeState(state.crossrange, 'length', true, constants);
    state.crossrange             = calculateCrossrange(state, constants);
    state.deadbandWidth          = calculateDeadbandWidth(state.velocity, constants);
end
