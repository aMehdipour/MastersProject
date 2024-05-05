% Predictor-corrector guidance function
function [commandedBankAngle, predictedTrajectory] = predictorCorrectorGuidance(parameters, constants, state, currentEnergy, finalEnergy, time)

    cntIterations = 0;

    bankAnglePertubation = deg2rad(1e-6);

    bankAngleHistory(cntIterations + 1) = abs(state.bankAngle);
    

    % Integrate equations of motion
    predictedTrajectory = integrateEquationsOfMotion(parameters,...
                                                     constants,...
                                                     state,...
                                                     currentEnergy,...
                                                     finalEnergy,...
                                                     bankAngleHistory(cntIterations + 1));

    % Evaluate terminal constraints
    rangeErrorHistory(cntIterations + 1) = evaluateTerminalConstraints(predictedTrajectory, parameters, constants);
    rangeErrorOld = 99999;

    % Iteratively adjust bank angle profile
    % while abs(rangeErrorHistory(cntIterations + 1)) > parameters.RANGE_ERROR_TOLERANCE / constants.EARTH_RADIUS_EQ...
    %    && abs(rangeErrorHistory(cntIterations + 1) - rangeErrorOld) > 1e-6
    while abs(rangeErrorHistory(cntIterations + 1) - rangeErrorOld) > 1e-9
        if cntIterations == 0
            % For the first pass, we use the Newton-Rhapson method to calculate the modified
            % bank angle. This is because we do not yet have the k-1 index necessary to use
            % the secant method.
            bankAngleHistory(cntIterations + 2) = bankAngleHistory(cntIterations + 1) + bankAnglePertubation;

            % Re-integrate equations of motion

            predictedTrajectory = integrateEquationsOfMotion(parameters,...
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
            predictedTrajectory = integrateEquationsOfMotion(parameters,...
                                                            constants,...
                                                            state,...
                                                            currentEnergy,...
                                                            finalEnergy,...
                                                            bankAngleHistory(cntIterations + 1));

            % Re-evaluate terminal constraint
            rangeErrorHistory(cntIterations + 1) = evaluateTerminalConstraints(predictedTrajectory, parameters, constants);
            rangeErrorOld = rangeErrorHistory(cntIterations + 1);
        end
    end

    for i = 1:length(predictedTrajectory.heading)
        [~, terminalAzimuth] = distance(predictedTrajectory.latitude(i), predictedTrajectory.longitude(i), deg2rad(constants.TARGET_LAT), deg2rad(constants.TARGET_LON));
        predictedTrajectory.headingError(i) = wrapToPi(terminalAzimuth - predictedTrajectory.heading(i));
    end

    crossrange = state.crossrange;
    energyStep = (finalEnergy - currentEnergy) / 20000;
    % Update the bank angle sign based on the lateral logic
    % bankAngleSign = sign(bankAngleHistory(end));

    if time == 0
        bankAngleSign = -sign(state.crossrange);
    end

    % for i = 1:length(predictedTrajectory.velocity)
    %     deadbandWidth = calculateDeadbandWidth(predictedTrajectory.velocity(i), constants);
    % 
    %     crossrangeRate = (cos(predictedTrajectory.flightPathAngle(i)) * sin(predictedTrajectory.headingError(i))) / (predictedTrajectory.normGeocentricDistance(i) * predictedTrajectory.drag(i));
    %     crossrange = crossrange + crossrangeRate * energyStep;
    % 
    %     if abs(crossrange) > deadbandWidth
    %         bankAngleSign = -bankAngleSign;
    %     end
    %     crossrangeHistory(i) = crossrange;
    %     deadbandWidthHistory(i) = deadbandWidth;
    % end

    % Return commanded bank angle
    commandedBankAngle = state.bankSign * bankAngleHistory(end);
end
