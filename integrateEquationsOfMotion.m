% Integrate equations of motion function
function predictedTrajectory = integrateEquationsOfMotion(parameters, constants, state, currentEnergy, finalEnergy, initialBankAngle)

    % Integration step size
    energyStep = (finalEnergy - currentEnergy) / 10000;
    
    % Initialize predicted trajectory
    energy = currentEnergy:energyStep:finalEnergy;
    rangeToGo = NaN(size(energy));
    normGeocentricDistance = NaN(size(energy));
    flightPathAngle = NaN(size(energy));
    liftHistory     = NaN(size(energy));
    dragHistory     = NaN(size(energy));
    velocityHistory = NaN(size(energy));
    altitudeHistory = NaN(size(energy));
    rhoHistory      = NaN(size(energy));
    trimClHistory   = NaN(size(energy));
    trimCdHistory   = NaN(size(energy));

    % Find the bank angle profile (Eq. 24 in Low Lifting Entry Guidance)
    bankAngleProfile = createBankAngleProfile(initialBankAngle, constants.FINAL_BANK_ANGLE, energy);
    
    % Set initial conditions
    rangeToGo(1) = state.rangeToGo / constants.EARTH_RADIUS_EQ;
    normGeocentricDistance(1) = state.r;
    flightPathAngle(1) = state.flightPathAngle;
    
    % Integrate equations of motion
    for i = 1:length(energy)-1
        % Compute velocity
        velocity = sqrt(2 * (1 / normGeocentricDistance(i) - energy(i)));

        if parameters.debug && ~isreal(velocity)
            warning("Imaginary velocity detected: \n Altitude: %.2d\n Energy: %.2d\n Range: %.2d\n Lift: %.2d\n Drag: %.2d", altitudeUnnormalized, energy(i), rangeToGo(i), lift, drag)
            break;
        end

        % Change states back to their dimensional forms
        velocityUnnormalized = velocity * constants.VELOCITY_SCALE;
        altitudeUnnormalized = (normGeocentricDistance(i) * constants.EARTH_RADIUS_EQ) - constants.EARTH_RADIUS_EQ;
        
        % Compute non-dimensionalized lift and drag forces
        [~, trimCD, trimCL] = calculateTrimAero(parameters, constants, velocityUnnormalized, altitudeUnnormalized);
        [rho, ~, ~] = atmosphereModel(altitudeUnnormalized);

        lift = 0.5 * rho * trimCL * velocityUnnormalized^2 / constants.VEHICLE_WEIGHT;
        drag = 0.5 * rho * trimCD * velocityUnnormalized^2 / constants.VEHICLE_WEIGHT;

        % Equations of motion
        rangeToGoRate = -cos(flightPathAngle(i)) / (normGeocentricDistance(i) * drag);
        normGeocentricDistanceRate = sin(flightPathAngle(i)) / drag;
        flightPathAngleRate = (lift * cos(bankAngleProfile(i)) + (velocity^2 - 1/normGeocentricDistance(i)) * (cos(flightPathAngle(i)) / normGeocentricDistance(i))) / (drag * velocity^2);
        
        % Update predicted trajectory
        rangeToGo(i+1) = rangeToGo(i) + rangeToGoRate * energyStep;
        normGeocentricDistance(i+1) = normGeocentricDistance(i) + normGeocentricDistanceRate * energyStep;
        flightPathAngle(i+1) = flightPathAngle(i) + flightPathAngleRate * energyStep;

        if parameters.debug && rangeToGo(i+1) <= 0.0
            warning("Range is less than zero: \n Altitude: %.2d\n Energy: %.2d\n Range: %.2d\n Lift: %.2d\n Drag: %.2d", altitudeUnnormalized, energy(i), rangeToGo(i), lift, drag)
            % break;
        end

        % Debug
        liftHistory(i) = lift;
        dragHistory(i) = drag;
        velocityHistory(i) = velocityUnnormalized;
        altitudeHistory(i) = altitudeUnnormalized;
        rhoHistory(i) = rho;
        trimClHistory(i) = trimCL;
        trimCdHistory(i) = trimCD;
    end
    
    % Debug plots
    if parameters.debug
        figure('Name','Normalized Lift And Drag vs Range')
        hold on
        plot(rangeToGo, dragHistory)
        plot(rangeToGo, liftHistory)
        plot(rangeToGo, sqrt(liftHistory.^2 + dragHistory.^2))
        xlabel('Normalized Range to Go')
        ylabel('Lift and Drag Accelerations (g)')
        legend('Lift', 'Drag', 'Total Load')

        figure('Name', 'Lift And Drag Coefficients')
        hold on
        plot(rangeToGo, trimClHistory)
        plot(rangeToGo, trimCdHistory)
        xlabel('Normalized Range to Go')
        ylabel('Trim Lift and Drag Coefficients')
        legend('$C_L$', '$C_D$')

        figure('Name', 'Altitude vs Range To Go')
        plot(rangeToGo, altitudeHistory)
        xlabel('Normalized Range to Go')
        ylabel('Altitude (m)')

        figure('Name', 'Velocity vs Range To Go')
        plot(rangeToGo, velocityHistory)
        xlabel('Normalized Range to Go')
        ylabel('Velocity (m/s)')

        figure('Name', 'FPA vs Range To Go')
        plot(rangeToGo, flightPathAngle)
        xlabel('Normalized Range to Go')
        ylabel('Flight Path Angle')

        figure('Name', 'Density History')
        plot(rangeToGo, rhoHistory)
        xlabel('Normalized Range to Go')
        ylabel('Density $(kg/m^3)$')

    end
    
    % Store predicted trajectory
    predictedTrajectory.energy = energy;
    predictedTrajectory.rangeToGo = rangeToGo;
    predictedTrajectory.normGeocentricDistance = normGeocentricDistance;
    predictedTrajectory.flightPathAngle = flightPathAngle;
    predictedTrajectory.bankAngleProfile = bankAngleProfile;
end
