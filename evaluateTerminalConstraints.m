% Evaluate terminal constraints function
function rangeError = evaluateTerminalConstraints(predictedTrajectory, parameters)
    % Extract the range-to-go at the final energy as determined
    % by the predictor step
    rangeToGoFinal = predictedTrajectory.range(end);
    
    % Calculate the range error (Eq. 20 in Entry Guidance: A Unified Method)
    rangeError =  rangeToGoFinal - parameters.RANGE_ERROR_TOLERANCE;
end

