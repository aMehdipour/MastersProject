% Integrate equations of motion function
function predictedTrajectory = integrateEquationsOfMotionFullURC(parameters, constants, state, currentEnergy, finalEnergy, initialAlpha, initialBeta)

    % Integration step size
    energyStep = (finalEnergy - currentEnergy) / 20000;

    % Initialize predicted trajectory
    energy                   = currentEnergy:energyStep:finalEnergy;
    rangeToGo                = NaN(size(energy));
    normGeocentricDistance   = NaN(size(energy));
    flightPathAngle          = NaN(size(energy));
    liftHistory              = NaN(size(energy));
    dragHistory              = NaN(size(energy));
    velocityHistory          = NaN(size(energy));
    altitudeHistory          = NaN(size(energy));
    rhoHistory               = NaN(size(energy));
    trimClHistory            = NaN(size(energy));
    trimCdHistory            = NaN(size(energy));
    heatingRateHistory       = NaN(size(energy));
    modifiedAlphaProfile     = NaN(size(energy)); % Changed to alpha profile
    crossrange               = NaN(size(energy));
    heading                  = NaN(size(energy));
    latitude                 = NaN(size(energy));
    longitude                = NaN(size(energy));
    deadband                 = NaN(size(energy));
    headingError             = NaN(size(energy));
    betaProfile              = NaN(size(energy));

    % Find the alpha profile (similar to Eq. 24 in Low Lifting Entry Guidance)
    alphaProfile = createAlphaProfile(initialAlpha, constants.FINAL_ALPHA, energy);

    % Set initial conditions
    rangeToGo(1)              = state.rangeToGo;
    normGeocentricDistance(1) = state.normGeocentricDistance;
    flightPathAngle(1)        = state.flightPathAngle;
    crossrange(1)             = state.crossrange;
    heading(1)                = state.heading;
    latitude(1)               = state.latitude;
    longitude(1)              = state.longitude;
    headingError(1)           = state.headingError;
    betaProfile(1)            = state.sideslipAngle;

    % Integrate equations of motion
    for i = 1:length(energy)-1
        % Compute velocity
        velocity = sqrt(2 * (1 / normGeocentricDistance(i) - energy(i)));

        % Change states back to their dimensional forms
        velocityUnnormalized = velocity * constants.VELOCITY_SCALE;
        altitudeUnnormalized = (normGeocentricDistance(i) * constants.EARTH_RADIUS_EQ) - constants.EARTH_RADIUS_EQ;

        % Compute non-dimensionalized lift, drag, and side forces
        [trimCL, trimCD, trimCS] = calculateAeroCoefficients(parameters, alphaProfile(i), betaProfile(i));
        [rho, ~, ~] = atmosphereModel(altitudeUnnormalized);

        lift = 0.5 * rho * trimCL * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;
        drag = 0.5 * rho * trimCD * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;
        sideForce = 0.5 * rho * trimCS * velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;

        % Equations of motion
        earthRotationRate                   = normalizeState(constants.EARTH_ROTATION_RATE, 'time', true, constants);
        flightPathAngleRate                 = (lift * cos(constants.INITIAL_BANK_ANGLE) - sideForce * sin(constants.INITIAL_BANK_ANGLE) + (velocity^2 - 1 / normGeocentricDistance(i)) * cos(flightPathAngle(i)) / normGeocentricDistance(i) + 2 * earthRotationRate * velocity * cos(latitude(i)) * sin(heading(i)) + earthRotationRate^2 * normGeocentricDistance(i) * cos(latitude(i)) * (cos(flightPathAngle(i)) * cos(latitude(i)) + sin(flightPathAngle(i)) * cos(heading(i)) * sin(latitude(i)))) / velocity;
        normGeocentricDistanceRate          = velocity * sin(flightPathAngle(i));
        rangeToGoRate                       = -velocity * cos(flightPathAngle(i)) / normGeocentricDistance(i);
        longitudeRate                       = (velocity * cos(flightPathAngle(i)) * sin(heading(i))) / (normGeocentricDistance(i) * cos(longitude(i)));
        latitudeRate                        = (velocity * cos(flightPathAngle(i)) * cos(heading(i))) / normGeocentricDistance(i);
        headingRate                         = 1 / velocity * ((lift * sin(constants.INITIAL_BANK_ANGLE)) / cos(flightPathAngle(i)) + (sideForce * sin(constants.INITIAL_BANK_ANGLE)) + (velocity^2 * cos(flightPathAngle(i)) * sin(heading(i)) * tan(latitude(i))) / normGeocentricDistance(i) - 2 * earthRotationRate * velocity * (tan(flightPathAngle(i)) * cos(heading(i)) * cos(latitude(i)) - sin(latitude(i))) + (earthRotationRate^2 * normGeocentricDistance(i) * sin(heading(i)) * sin(latitude(i)) * cos(latitude(i))) / cos(flightPathAngle(i)) );
        [~, terminalAzimuth] = distance('gc',latitude(i), longitude(i), deg2rad(constants.TARGET_LAT), deg2rad(constants.TARGET_LON), 'radians');
        headingError(i+1) = wrapToPi(terminalAzimuth - heading(i));



        % Update predicted trajectory
        rangeToGo(i+1) = rangeToGo(i) + rangeToGoRate / (drag * velocity) * energyStep;
        normGeocentricDistance(i+1) = normGeocentricDistance(i) + normGeocentricDistanceRate / (drag * velocity) * energyStep;
        flightPathAngle(i+1) = flightPathAngle(i) + flightPathAngleRate / (drag * velocity) * energyStep;
        heading(i+1) = heading(i) + headingRate / (drag * velocity) * energyStep;
        latitude(i+1) = latitude(i) + latitudeRate / (drag * velocity) * energyStep;
        longitude(i+1) = longitude(i) + longitudeRate / (drag * velocity) * energyStep;

        % Calculate the lateral logic
        gamma = constants.VEHICLE_MASS / (0.5 * rho * velocityUnnormalized^2 * constants.FRONTAL_AREA);
        maxSideForce = max(max(abs(parameters.CStable)));
        headingGain = abs((maxSideForce / gamma) / (0.5 * maxSideForce));
        headingRateGain = 2 * 1 * constants.VEHICLE_MASS * sqrt(headingGain / constants.VEHICLE_MASS);

        sideForceCommand = gamma * (headingGain * (terminalAzimuth - heading(i)) + headingRateGain * (0 - headingRate));

        betaProfile(i+1) = findSideslipAngle(parameters,alphaProfile(i),sideForceCommand);

        % Debug
        liftHistory(i) = lift;
        dragHistory(i) = drag;
        velocityHistory(i) = velocityUnnormalized;
        altitudeHistory(i) = altitudeUnnormalized;
        rhoHistory(i) = rho;
        trimClHistory(i) = trimCL;
        trimCdHistory(i) = trimCD;
        heatingRate = constants.K_Q * sqrt(rho) * velocity^3.15 * 1e-4;
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
    predictedTrajectory.alphaProfile = alphaProfile;
    predictedTrajectory.heading = heading;
    predictedTrajectory.latitude = latitude;
    predictedTrajectory.longitude = longitude;
    predictedTrajectory.velocity = velocityHistory;
    predictedTrajectory.crossrange = crossrange;
    predictedTrajectory.drag = dragHistory;
    predictedTrajectory.beta = betaProfile;

end
