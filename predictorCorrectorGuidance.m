% Predictor-corrector guidance function
function [commandedBankAngle, predictedTrajectory] = predictorCorrectorGuidance(parameters, constants, state, currentEnergy, finalEnergy, time)

    cntIterations = 0;

    bankAnglePertubation = deg2rad(1e-3);

    bankAngleHistory(cntIterations + 1) = abs(state.bankAngle);
    

    % Integrate equations of motion
    predictedTrajectory = integrateEquationsOfMotionFull(parameters,...
                                                     constants,...
                                                     state,...
                                                     currentEnergy,...
                                                     finalEnergy,...
                                                     bankAngleHistory(cntIterations + 1));

    % Evaluate terminal constraints
    rangeErrorHistory(cntIterations + 1) = evaluateTerminalConstraints(predictedTrajectory, parameters, constants);
    rangeErrorOld = 99999;

    % Iteratively adjust bank angle profile
    while abs(rangeErrorHistory(cntIterations + 1) - rangeErrorOld) > 1e-3...
       && cntIterations < 20
        if cntIterations == 0
            % For the first pass, we use the Newton-Rhapson method to calculate the modified
            % bank angle. This is because we do not yet have the k-1 index necessary to use
            % the secant method.
            bankAngleHistory(cntIterations + 2) = bankAngleHistory(cntIterations + 1) + bankAnglePertubation;

            % Re-integrate equations of motion

            predictedTrajectory = integrateEquationsOfMotionFull(parameters,...
                                                             constants,...
                                                             state,...
                                                             currentEnergy,...
                                                             finalEnergy,...
                                                             bankAngleHistory(cntIterations + 2));

            % Re-evaluate terminal constraint
            rangeErrorHistory(cntIterations + 2) = evaluateTerminalConstraints(predictedTrajectory, parameters, constants);

            % Adjust bank angle profile based on terminal constraint
            bankAngleHistory(cntIterations + 3) = modifyBankAngleProfile(rangeErrorHistory, bankAngleHistory, cntIterations);

            cntIterations = cntIterations + 1;
        else

            % Adjust bank angle profile based on terminal constraint
            bankAngleHistory(cntIterations + 2) = modifyBankAngleProfile(rangeErrorHistory, bankAngleHistory, cntIterations);

            cntIterations = cntIterations + 1;

            % Re-integrate equations of motion
            predictedTrajectory = integrateEquationsOfMotionFull(parameters,...
                                                            constants,...
                                                            state,...
                                                            currentEnergy,...
                                                            finalEnergy,...
                                                            bankAngleHistory(cntIterations + 1));

            % Re-evaluate terminal constraint
            rangeErrorHistory(cntIterations + 1) = evaluateTerminalConstraints(predictedTrajectory, parameters, constants);
            rangeErrorOld = rangeErrorHistory(cntIterations);
        end
    end

    % Return commanded bank angle
    idx = find(abs(rangeErrorHistory) == min(abs(rangeErrorHistory)));
    commandedBankAngle = state.bankSign * abs(bankAngleHistory(end));
end
