function [state, derivatives] = updateStateEnergy(state, commandedBankAngle, constants, parameters, energyStep)
    velocity               = state.velocity;
    flightPathAngle        = state.flightPathAngle;
    bankAngle              = state.bankAngle;
    latitude               = state.latitude;
    longitude              = state.longitude;
    heading                = state.heading;
    normGeocentricDistance = state.normGeocentricDistance;
    rangeToGo              = state.rangeToGo;
    velocityUnnormalized = state.velocityUnnormalized;
    altitudeUnnormalized = state.altitudeUnnormalized;

    % Compute non-dimensionalized lift and drag forces
    [rho, ~, ~] = atmosphereModel(altitudeUnnormalized);
    [~, trimCD, trimCL] = calculateTrimAero(parameters, constants, velocityUnnormalized, altitudeUnnormalized);
    lift = 0.5 * rho * trimCL * velocityUnnormalized^2 / constants.VEHICLE_WEIGHT;
    drag = 0.5 * rho * trimCD * velocityUnnormalized^2 / constants.VEHICLE_WEIGHT;

    % Equations of motion
    % derivatives.velocityDot               = -drag - sin(flightPathAngle) / normGeocentricDistance^2;
    derivatives.normGeocentricDistanceDot = sin(flightPathAngle) / drag;
    derivatives.flightPathAngleDot        = (lift * cos(commandedBankAngle) + (velocity^2 - 1 / normGeocentricDistance) * cos(flightPathAngle) / normGeocentricDistance) / (drag * velocity^2);
    derivatives.rangeDot                  = -cos(flightPathAngle) / (normGeocentricDistance * drag);
    derivatives.lonDot                    = (cos(flightPathAngle) * sin(heading)) / (normGeocentricDistance * cos(longitude) * drag);
    derivatives.latDot                    = (cos(flightPathAngle) * cos(heading)) / (normGeocentricDistance * drag);
    derivatives.headingDot                = ((lift * sin(bankAngle)) / cos(flightPathAngle) + (velocity^2) / normGeocentricDistance * cos(flightPathAngle) * sin(heading) * tan(latitude)) / (velocity * drag);
    % derivatives.bankAngleDot              = clampBankAngleRate(bankAngle, commandedBankAngle, constants, energyStep);

    % Update state variables
    % state.velocity               = velocity + derivatives.velocityDot * energyStep;
    state.normGeocentricDistance = normGeocentricDistance + derivatives.normGeocentricDistanceDot * energyStep;
    state.flightPathAngle        = flightPathAngle + derivatives.flightPathAngleDot * energyStep;
    state.rangeToGo              = rangeToGo + derivatives.rangeDot * energyStep;
    state.longitude              = longitude + derivatives.lonDot * energyStep;
    state.latitude               = latitude + derivatives.latDot * energyStep;
    state.heading                = heading + derivatives.headingDot * energyStep;
    % state.bankAngle              = bankAngle + derivatives.bankAngleDot;
    state.bankAngle              = commandedBankAngle;
    state.altitudeUnnormalized   = normalizeState(state.normGeocentricDistance,'length',false,constants) - constants.EARTH_RADIUS_EQ;
    state.velocityUnnormalized   = normalizeState(state.velocity,'velocity',false,constants);
    state.lift                   = lift;
    state.drag                   = drag;
    state.loadFactor             = sqrt(lift^2 + drag^2);
    state.heatingRate            = constants.K_Q * sqrt(rho) * velocity^3.15 * 1e-4;
    state.altitudeUnnormalized   = normalizeState(state.normGeocentricDistance,'length',false,constants) - constants.EARTH_RADIUS_EQ;
    state.velocityUnnormalized   = normalizeState(state.velocity,'velocity',false,constants);
    state.rangeToGoUnnormalized  = normalizeState(state.rangeToGo, 'length', false, constants);
end
