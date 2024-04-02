% Predictor-corrector guidance function
function [commandedBankAngle, predictedTrajectory] = predictorCorrectorGuidance(parameters, constants, r, currentEnergy, finalEnergy, currentRange, bankAngle)

    cntIterations = 0;

    bankAngleHistory(cntIterations + 1) = bankAngle;

    % Integrate equations of motion
    predictedTrajectory = integrateEquationsOfMotion(parameters,...
                                                     constants,...
                                                     r,...
                                                     currentEnergy,...
                                                     finalEnergy,...
                                                     currentRange,...
                                                     bankAngleHistory(cntIterations + 1));

    % Evaluate terminal constraints
    rangeErrorHistory(cntIterations + 1) = evaluateTerminalConstraints(predictedTrajectory, parameters);

    % Iteratively adjust bank angle profile
    while abs(rangeErrorHistory(cntIterations + 1)) > 1e-3
        if cntIterations == 0
            % For the first pass, we use the Newton-Rhapson method to calculate the modified
            % bank angle. This is because we do not yet have the k-1 index necessary to use
            % the secant method.
            bankAngleHistory(cntIterations + 2) = bankAngle + epsilon;

            % Re-integrate equations of motion

            predictedTrajectory = integrateEquationsOfMotion(parameters,...
                                                             constants,...
                                                             r,...
                                                             currentEnergy,...
                                                             finalEnergy,...
                                                             currentRange,...
                                                             bankAngleHistory(cntIterations + 2));

            % Re-evaluate terminal constraint
            rangeErrorHistory = evaluateTerminalConstraints(predictedTrajectory, parameters);

            % Adjust bank angle profile based on terminal constraint
            bankAngleHistory(cntIterations + 3) = modifyBankAngleProfile(rangeErrorHistory, bankAngleHistory, cntIterations);

            cntIterations = cntIterations + 1;

            continue
        end

        % Adjust bank angle profile based on terminal constraint
        bankAngleHistory(cntIterations + 2) = modifyBankAngleProfile(rangeErrorHistory, bankAngleHistory, cntIterations);

        cntIterations = cntIterations + 1;

        % Re-integrate equations of motion
        predictedTrajectory = integrateEquationsOfMotion(parameters,...
                                                         constants,...
                                                         r,...
                                                         currentEnergy,...
                                                         finalEnergy,...
                                                         currentRange,...
                                                         bankAngleHistory(cntIterations + 1));

        % Re-evaluate terminal constraint
        rangeErrorHistory = evaluateTerminalConstraints(predictedTrajectory, parameters);
    end

    % Return commanded bank angle
    commandedBankAngle = initialBankAngle;
end
