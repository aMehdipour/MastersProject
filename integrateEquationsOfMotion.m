% Integrate equations of motion function
function predictedTrajectory = integrateEquationsOfMotion(parameters, constants, r, currentEnergy, finalEnergy, currentRange, initialBankAngle)

    % Integration step size
    energyStep = (finalEnergy - currentEnergy) / 1000;
    
    % Initialize predicted trajectory
    energy = currentEnergy:energyStep:finalEnergy;
    range = zeros(size(energy));
    altitude = zeros(size(energy));
    flightPathAngle = zeros(size(energy));

    % Find the bank angle profile (Eq. 24 in Low Lifting Entry Guidance)
    bankAngleProfile = createBankAngleProfile(initialBankAngle, constants.FINAL_BANK_ANGLE, energy);
    
    % Set initial conditions
    range(1) = currentRange;
    altitude(1) = r;
    flightPathAngle(1) = 0;
    
    % Integrate equations of motion
    for i = 1:length(energy)-1
        % Compute velocity
        velocity = sqrt(2 * (1 / altitude(i) - energy(i)));
        velocityUnnormalized = velocity * constants.VELOCITY_SCALE;
        
        % Compute non-dimensionalized lift and drag forces
        [~, trimCD, trimCL] = calculateTrimAero(parameters, constants, velocityUnnormalized, altitude(i));
        [rho, ~, ~] = atmosphereModel(H0);

        lift = 0.5 * rho * trimCL * velocity^2 / constants.VEHICLE_WEIGHT;
        drag = 0.5 * rho * trimCD * velocity^2 / constants.VEHICLE_WEIGHT;


        % Equations of motion
        rangeRate = -cos(flightPathAngle(i)) / (altitude(i) * drag);
        altitudeRate = sin(flightPathAngle(i)) / drag;
        flightPathAngleRate = (lift * cos(bankAngleProfile(i)) + (velocity^2 - 1/altitude(i)) * (cos(flightPathAngle(i)) / altitude(i))) / (drag * velocity^2);
        
        % Update predicted trajectory
        range(i+1) = range(i) + rangeRate * energyStep;
        altitude(i+1) = altitude(i) + altitudeRate * energyStep;
        flightPathAngle(i+1) = flightPathAngle(i) + flightPathAngleRate * energyStep;
    end
    
    % Store predicted trajectory
    predictedTrajectory.energy = energy;
    predictedTrajectory.range = range;
    predictedTrajectory.altitude = altitude;
    predictedTrajectory.flightPathAngle = flightPathAngle;
    predictedTrajectory.bankAngleProfile = bankAngleProfile;
end
