% Predictor-corrector guidance function
function [commandedAlpha, commandedBeta, predictedTrajectory] = predictorCorrectorGuidanceURC(parameters, constants, state, currentEnergy, finalEnergy, time)

    cntIterations = 0;

    alphaPerturbation = deg2rad(1);

    alphaHistory(cntIterations + 1) = abs(state.angleOfAttack);
    betaHistory(cntIterations + 1) = state.sideslipAngle;
    
    % Integrate equations of motion
    predictedTrajectory = integrateEquationsOfMotionFullURC(parameters,...
                                                     constants,...
                                                     state,...
                                                     currentEnergy,...
                                                     finalEnergy,...
                                                     alphaHistory(cntIterations + 1));

    % Evaluate terminal constraints
    rangeErrorHistory(cntIterations + 1) = evaluateTerminalConstraints(predictedTrajectory, parameters, constants);
    rangeErrorOld = 99999;

    % Iteratively adjust alpha profile
    while abs(rangeErrorHistory(cntIterations + 1) - rangeErrorOld) > 1e-3...
       && cntIterations < 20
        if cntIterations == 0
            % For the first pass, we use the Newton-Rhapson method to calculate the modified
            % alpha. This is because we do not yet have the k-1 index necessary to use
            % the secant method.
            alphaHistory(cntIterations + 2) = alphaHistory(cntIterations + 1) + alphaPerturbation;

            % Re-integrate equations of motion
            predictedTrajectory = integrateEquationsOfMotionFullURC(parameters,...
                                                             constants,...
                                                             state,...
                                                             currentEnergy,...
                                                             finalEnergy,...
                                                             alphaHistory(cntIterations + 2));

            % Re-evaluate terminal constraint
            rangeErrorHistory(cntIterations + 2) = evaluateTerminalConstraints(predictedTrajectory, parameters, constants);

            % Adjust alpha profile based on terminal constraint
            alphaHistory(cntIterations + 3) = modifyAlphaProfile(rangeErrorHistory, alphaHistory, cntIterations);

            cntIterations = cntIterations + 1;
        else
            % Adjust alpha profile based on terminal constraint
            alphaHistory(cntIterations + 2) = modifyAlphaProfile(rangeErrorHistory, alphaHistory, cntIterations);

            cntIterations = cntIterations + 1;

            % Re-integrate equations of motion
            predictedTrajectory = integrateEquationsOfMotionFullURC(parameters,...
                                                            constants,...
                                                            state,...
                                                            currentEnergy,...
                                                            finalEnergy,...
                                                            alphaHistory(cntIterations + 1));

            % Re-evaluate terminal constraint
            rangeErrorHistory(cntIterations + 1) = evaluateTerminalConstraints(predictedTrajectory, parameters, constants);
            rangeErrorOld = rangeErrorHistory(cntIterations);
        end
    end

    % Return commanded alpha and beta
    idx = find(abs(rangeErrorHistory) == min(abs(rangeErrorHistory)));
    commandedAlpha = alphaHistory(idx);
    commandedBeta = predictedTrajectory.beta(2);
end
