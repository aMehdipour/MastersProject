% Integrate equations of motion function
function predictedTrajectory = integrateEquationsOfMotionFull(parameters, constants, state, currentEnergy, finalEnergy, initialBankAngle)

    % Integration step size
    energyStep = (finalEnergy - currentEnergy) / 20000;

    % Initialize predicted trajectory
    energy                   = currentEnergy:energyStep:finalEnergy;
    rangeToGo                = NaN(size(energy));
    normGeocentricDistance   = NaN(size(energy));
    flightPathAngle          = NaN(size(energy));
    flightPathAngleTest      = NaN(size(energy));
    liftHistory              = NaN(size(energy));
    dragHistory              = NaN(size(energy));
    velocityHistory          = NaN(size(energy));
    altitudeHistory          = NaN(size(energy));
    rhoHistory               = NaN(size(energy));
    trimClHistory            = NaN(size(energy));
    trimCdHistory            = NaN(size(energy));
    heatingRateHistory       = NaN(size(energy));
    modifiedBankAngleProfile = NaN(size(energy));
    crossrange               = NaN(size(energy));
    heading                  = NaN(size(energy));
    latitude                 = NaN(size(energy));
    longitude                = NaN(size(energy));
    deadband                 = NaN(size(energy));
    headingError             = NaN(size(energy));

    % Find the bank angle profile (Eq. 24 in Low Lifting Entry Guidance)
    bankAngleProfile = createBankAngleProfile(initialBankAngle, constants.FINAL_BANK_ANGLE, energy);

    % Set initial conditions
    rangeToGo(1)              = state.rangeToGo;
    normGeocentricDistance(1) = state.normGeocentricDistance;
    flightPathAngle(1)        = state.flightPathAngle;
    crossrange(1)             = state.crossrange;
    heading(1)                = state.heading;
    latitude(1)               = state.latitude;
    longitude(1)              = state.longitude;
    headingError(1)           = state.headingError;

    bankAngleSign = state.bankSign;

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

        lift = 0.5 * rho * trimCL * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;
        drag = 0.5 * rho * trimCD * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;

        bankAngleCommandTest = evaluatePathConstraints(altitudeUnnormalized, velocityUnnormalized, rho, lift, drag, velocity, energyStep, flightPathAngle(i), bankAngleProfile(i), constants);

        % if parameters.debug && ~isreal(bankAngleCommand)
        %     warning("Imaginary bank angle detected: \n Altitude: %.2d\n Energy: %.2d\n Range: %.2d\n Lift: %.2d\n Drag: %.2d", altitudeUnnormalized, energy(i), rangeToGo(i), lift, drag)
        %     break;
        % end

        % if isreal(bankAngleCommand)
            modifiedBankAngleProfile(i) = bankAngleCommandTest;
        % end

        % Equations of motion
        % rangeToGoRate = -cos(flightPathAngle(i)) / (normGeocentricDistance(i) * drag);
        % normGeocentricDistanceRate = sin(flightPathAngle(i)) / drag;
        % flightPathAngleRate = (lift * cos(bankAngleProfile(i)) + (velocity^2 - 1/normGeocentricDistance(i)) * (cos(flightPathAngle(i)) / normGeocentricDistance(i))) / (drag * velocity^2);
        % headingRate = (lift * sin(bankAngleProfile(i))) / (drag * velocity^2 * cos(flightPathAngle(i))) + cos(flightPathAngle(i)) * sin(heading(i)) * tan(latitude(i)) / (normGeocentricDistance(i) * drag);
        % latitudeRate = (cos(flightPathAngle(i)) * cos(heading(i))) / (normGeocentricDistance(i) * drag);
        % longitudeRate = (cos(flightPathAngle(i)) * sin(heading(i))) / (normGeocentricDistance(i) * drag * cos(latitude(i)));
        earthRotationRate                   = normalizeState(constants.EARTH_ROTATION_RATE, 'time', true, constants);
        flightPathAngleRate                 = (lift * cos(bankAngleProfile(i)) + (velocity^2 - 1 / normGeocentricDistance(i)) * cos(flightPathAngle(i)) / normGeocentricDistance(i) + 2 * earthRotationRate * velocity * cos(latitude(i)) * sin(heading(i)) + earthRotationRate^2 * normGeocentricDistance(i) * cos(latitude(i)) * (cos(flightPathAngle(i)) * cos(latitude(i)) + sin(flightPathAngle(i)) * cos(heading(i)) * sin(latitude(i)))) / velocity;
        normGeocentricDistanceRate          = velocity * sin(flightPathAngle(i));
        rangeToGoRate                       = -velocity * cos(flightPathAngle(i)) / normGeocentricDistance(i);
        longitudeRate                       = (velocity * cos(flightPathAngle(i)) * sin(heading(i))) / (normGeocentricDistance(i) * cos(longitude(i)));
        latitudeRate                        = (velocity * cos(flightPathAngle(i)) * cos(heading(i))) / normGeocentricDistance(i);
        headingRate                         = 1 / velocity * ((lift * sin(bankAngleProfile(i))) / cos(flightPathAngle(i)) + (velocity^2 * cos(flightPathAngle(i)) * sin(heading(i)) * tan(latitude(i))) / normGeocentricDistance(i) - 2 * earthRotationRate * velocity * (tan(flightPathAngle(i)) * cos(heading(i)) * cos(latitude(i)) - sin(latitude(i))) + (earthRotationRate^2 * normGeocentricDistance(i) * sin(heading(i)) * sin(latitude(i)) * cos(latitude(i))) / cos(flightPathAngle(i)) );
        [~, terminalAzimuth] = distance('gc',latitude(i), longitude(i), deg2rad(constants.TARGET_LAT), deg2rad(constants.TARGET_LON), 'radians');
        headingError(i+1) = wrapToPi(terminalAzimuth - heading(i));

        % Update flight path angle rate using modified bank angle command
        flightPathAngleRateTest = (lift * cos(bankAngleCommandTest) + (velocity^2 - 1/normGeocentricDistance(i)) * (cos(flightPathAngle(i)) / normGeocentricDistance(i))) / (drag * velocity^2);

        % Update predicted trajectory
        stateTest.rangeToGo = rangeToGo(i);
        stateTest.heading = heading(i);
        stateTest.latitude = latitude(i);
        stateTest.longitude = longitude(i);
        rangeToGo(i+1) = rangeToGo(i) + rangeToGoRate / (drag * velocity) * energyStep;
        normGeocentricDistance(i+1) = normGeocentricDistance(i) + normGeocentricDistanceRate / (drag * velocity) * energyStep;
        flightPathAngle(i+1) = flightPathAngle(i) + flightPathAngleRate / (drag * velocity) * energyStep;
        flightPathAngleTest(i+1) = flightPathAngleTest(i) + flightPathAngleRateTest / (drag * velocity) * energyStep;
        crossrange(i+1) = calculateCrossrange(stateTest, constants);
        heading(i+1) = heading(i) + headingRate / (drag * velocity) * energyStep;
        latitude(i+1) = latitude(i) + latitudeRate / (drag * velocity) * energyStep;
        longitude(i+1) = longitude(i) + longitudeRate / (drag * velocity) * energyStep;

        deadband(i) = calculateDeadbandWidth(velocity,constants);

        if crossrange(i) > calculateDeadbandWidth(velocity,constants)
            bankAngleSign = -1;
        elseif crossrange(i) < -calculateDeadbandWidth(velocity,constants)
            bankAngleSign = 1;
        end

        bankAngleProfile(i+1) = bankAngleSign * bankAngleProfile(i+1);

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
        heatingRate = constants.K_Q * sqrt(rho) * velocity^3.15 * 1e-4;
        if heatingRate > 300
            % disp("Just about everywhere, it's gonna be hot")
        end
        heatingRateHistory(i) = heatingRate;
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

        figure('Name', 'Heating Rate History')
        plot(rangeToGo, heatingRateHistory)
        xlabel('Normalized Range to Go')
        ylabel('Heating Rate $(W/cm^2)$')

        figure('Name', 'Velocity vs Altitude')
        plot(velocityHistory, altitudeHistory./1e3)
        xlabel('Velocity m/s')
        ylabel('Altitude (km)')

        figure('Name', 'Latitude vs Range To Go')
        plot(rangeToGo, latitude)
        xlabel('Normalized Range to Go')
        ylabel('Latitude')

        figure('Name', 'Longitude vs Range To Go')
        plot(rangeToGo, longitude)
        xlabel('Normalized Range to Go')
        ylabel('Longitude')

        figure('Name', 'Heading vs Range To Go')
        plot(rangeToGo, heading)
        xlabel('Normalized Range to Go')
        ylabel('Heading')

        figure('Name', 'Crossrange vs Range To Go')
        plot(rangeToGo, crossrange)
        hold on
        plot(rangeToGo, deadband, '--r')
        plot(rangeToGo, -deadband, '--r')
        xlabel('Normalized Range to Go')
        ylabel('Crossrange')


        figure('Name', 'Bank Angle Profile Comparison')
        hold on
        plot(rangeToGo, rad2deg(bankAngleProfile))
        plot(rangeToGo, rad2deg(modifiedBankAngleProfile), '--')
        xlabel('Normalized Range To Go')
        ylabel('Bank Angle (deg)')

    end

    % Store predicted trajectory
    predictedTrajectory.energy = energy;
    predictedTrajectory.rangeToGo = rangeToGo;
    predictedTrajectory.normGeocentricDistance = normGeocentricDistance;
    predictedTrajectory.flightPathAngle = flightPathAngle;
    % predictedTrajectory.bankAngleProfile = modifiedBankAngleProfile;
    predictedTrajectory.bankAngleProfile = bankAngleProfile;
    predictedTrajectory.heading = heading;
    predictedTrajectory.latitude = latitude;
    predictedTrajectory.longitude = longitude;
    predictedTrajectory.velocity = velocityHistory;
    predictedTrajectory.crossrange = crossrange;
    predictedTrajectory.drag = dragHistory;

end
