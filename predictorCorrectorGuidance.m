% Predictor-corrector guidance function
function [commandedBankAngle, predictedTrajectory] = predictorCorrectorGuidance(parameters, constants, r, currentEnergy, finalEnergy, currentRange, bankAngle)
    
    cntIterations = 1;

    bankAngleHistory(cntIterations) = bankAngle;

    % Integrate equations of motion
    predictedTrajectory = integrateEquationsOfMotion(parameters, constants, r, currentEnergy, finalEnergy, currentRange, bankAngle);
    
    % Evaluate terminal constraints
    rangeError = evaluateTerminalConstraints(predictedTrajectory, parameters);
    
    % Iteratively adjust bank angle profile
    while abs(rangeError) > 1e-3
        if cntIterations == 1
            cntIterations = cntIterations + 1;

            bankAngle = bankAngle + epsilon;

            bankAngleHistory(cntIterations) = bankAngle;

            % Re-integrate equations of motion
            predictedTrajectory = integrateEquationsOfMotion(currentEnergy, finalEnergy, currentRange, bankAngleProfile, CL, CD);

            % Re-evaluate terminal constraint
            rangeError = evaluateTerminalConstraints(predictedTrajectory, parameters);

            % Adjust bank angle profile based on terminal constraint
            bankAngleHistory(cntIterations + 1) = modifyBankAngleProfile(rangeError, bankAngleHistory, cntIterations);

            continue
        end
        % Adjust bank angle profile based on terminal constraint
        bankAngleHistory(cntIterations + 1) = modifyBankAngleProfile(rangeError, bankAngleHistory, cntIterations);
        
        % Re-integrate equations of motion
        predictedTrajectory = integrateEquationsOfMotion(currentEnergy, finalEnergy, currentRange, bankAngleProfile, CL, CD);
        
        % Re-evaluate terminal constraint
        rangeError = evaluateTerminalConstraints(predictedTrajectory, parameters);
    end
    
    % Return commanded bank angle
    commandedBankAngle = initialBankAngle;
end
