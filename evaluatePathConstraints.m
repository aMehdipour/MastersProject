function bankAngleCommand = evaluatePathConstraints(altitudeUnnormalized, velocityUnnormalized, rho, lift, drag, velocity, energyStep, flightPathAngle, bankAngle, constants)
        % Calculate β_r using complex-step derivative approximation
        h = 1e-20; % Perturbation step size
        altitudeComplex = altitudeUnnormalized + 1i * h; % 1i here denotes the imaginary number
        rhoComplex = atmosphereModel(altitudeComplex);
        dRhoDr = imag(rhoComplex) / h;
        beta_r = -1000;%dRhoDr / rho * constants.EARTH_RADIUS_EQ;

        %% Load factor
        % Compute load factor and its energy-based derivative
        loadFactor = sqrt(lift^2 + drag^2);

        % Compute pseudocontrol requirement U_loadFactor
        % deltaEnergy = energyStep;
        deltaTime = 0.02;%16 / constants.TIME_SCALE;
        % deltaEnergy = drag * velocity * deltaTime;
        % A_loadFactor = -(2 * loadFactor * drag / velocity) / (drag * velocity);
        % B_loadFactor = (loadFactor * velocity * beta_r) / (drag * velocity);
        % U_loadFactor = (constants.MAX_LOAD_FACTOR - loadFactor - A_loadFactor * deltaEnergy) / (B_loadFactor * deltaEnergy);
        A_loadFactor = -(2 * loadFactor * drag / velocity);
        B_loadFactor = (loadFactor * velocity * beta_r);
        U_loadFactor = (constants.MAX_LOAD_FACTOR - loadFactor - A_loadFactor * deltaTime) / (B_loadFactor * deltaTime);
        U_loadFactor = (constants.MAX_LOAD_FACTOR - loadFactor * (1 - 2 * drag * deltaTime / velocity)) / (B_loadFactor * deltaTime);

        %% Heating rate
        % Calculate heating rate and its energy-based derivative
        heatingRate = constants.K_Q * sqrt(rho) * velocity^3.15;

        % Compute pseudocontrol requirement U_heatingRate
        % deltaEnergy = energyStep;
        deltaTime = 50 / constants.TIME_SCALE;
        deltaEnergy = drag * velocity * deltaTime;
        A_heatingRate = (3.15 * constants.K_Q * (rho * velocity)^1.575 * (-1.575 * drag * velocity)) / ...
               (drag * velocity);
        B_heatingRate = (3.15 * constants.K_Q * (rho * velocity)^1.575 * (0.5 * rho * velocity * beta_r)) / ...
               (drag * velocity);
        U_heatingRate = (constants.MAX_HEATING_RATE - heatingRate - A_heatingRate * deltaEnergy) / (B_heatingRate * deltaEnergy);
        % A_heatingRate = (3.15 * constants.K_Q * (rho * velocityUnnormalized)^1.575 * (-1.575 * drag * velocityUnnormalized));
        % B_heatingRate = (3.15 * constants.K_Q * (rho * velocityUnnormalized)^1.575 * (0.5 * rho * velocityUnnormalized * beta_r));
        % U_heatingRate = (constants.MAX_HEATING_RATE - heatingRate - A_heatingRate * deltaTime) / (B_heatingRate * deltaTime);

        % Calculate reference flight path angle and altitude rate
        if U_loadFactor > sin(flightPathAngle)
            test = 1;
        end
        
        sinGammaRef = max([sin(flightPathAngle), U_loadFactor, U_heatingRate]);
        refAltitudeRate = velocity * sinGammaRef;

        % Modify bank angle command using Eq. (41)
        k_0 = 100; % Control law gain
        altitudeRate = velocity * sin(flightPathAngle);
        bankAngleCommand = acos((lift * cos(bankAngle) - k_0 * (altitudeRate - refAltitudeRate)) / lift);

        if loadFactor > constants.MAX_LOAD_FACTOR
            % disp("THAT'S A BIG LOAD")
        end

        if ~isreal(bankAngleCommand)
             % warning('Imaginary bank angle command in evaluatePathConstraints!\n\tVelocity: %.2f\n\tAltitude: %.2f\n\tBank angle command:%.2f\n\tLoad Factor:%.2f\n\tHeating Rate:%.2f',...
             %                velocityUnnormalized, altitudeUnnormalized, bankAngleCommand, loadFactor, heatingRate)
        end
