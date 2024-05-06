% Evaluate terminal constraints function
function rangeError = evaluateTerminalConstraints(predictedTrajectory, parameters, constants)
    % Extract the range-to-go at the final energy as determined
    % by the predictor step
    rangeToGoFinal = predictedTrajectory.rangeToGo(end);
    % [rangeToGoFinal, ~] = distance('gc', predictedTrajectory.latitude(end), predictedTrajectory.longitude(end), deg2rad(constants.TARGET_LAT), deg2rad(constants.TARGET_LON));
    % rangeToGoFinal = rangeToGoFinal * sign(predictedTrajectory.rangeToGo(end));
    % Calculate the range error (Eq. 20 in Entry Guidance: A Unified Method)
    rangeError =  rangeToGoFinal - parameters.RANGE_ERROR_TOLERANCE / constants.EARTH_RADIUS_EQ;
end

