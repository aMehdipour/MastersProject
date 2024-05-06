clear, close all

% Constants
constants.GRAVITY_NOMINAL    = 9.81;           % Gravitational acceleration (m/s^2)
constants.EARTH_RADIUS_EQ    = 6378137;        % Earth's radius (m)
constants.INITIAL_ALTITUDE   = 122000;         % Initial altitude (m)
constants.INITIAL_VELOCITY   = 10100;          % Initial velocity (m/s)
constants.INITIAL_FPA        = -5.2 * pi/180;  % Initial flight path angle (rad)
constants.INITIAL_LAT        = -4.7;           % EI latitude (deg)
constants.INITIAL_LON        = -112;           % EI longitude (deg)
constants.TARGET_ALTITUDE    = 31000;          % Target final altitude (m)
constants.TARGET_VELOCITY    = 690;            % Target final velocity (m/s)
constants.TARGET_LAT         = 40;             % Final latitude (deg)
constants.TARGET_LON         = -112;           % Final longitude (deg)
constants.VEHICLE_MASS       = 75.7;           % Vehicle mass (kg)
constants.INITIAL_BANK_ANGLE = 0;
constants.FINAL_BANK_ANGLE   = deg2rad(70);
constants.REFERENCE_LENGTH   = 0.69; % reference length for moment calculation (m)
constants.VEHICLE_WEIGHT     = constants.VEHICLE_MASS * constants.GRAVITY_NOMINAL;
constants.FLATTENING_FACTOR  = 1 / 298.257223563;
constants.EARTH_RADIUS_POLE  = constants.EARTH_RADIUS_EQ * (1 - constants.FLATTENING_FACTOR);
constants.MAX_LOAD_FACTOR    = 3;             % Max load factor (g)
constants.K_Q                = 9.4369e-5 * (sqrt(constants.GRAVITY_NOMINAL * constants.EARTH_RADIUS_EQ))^3.15;
constants.BANK_ANGLE_RATE_LIMIT = deg2rad(10); % 10 degree bank angle rate limit
constants.MAX_HEATING_RATE = 2500000; % Maximum heating rate (w/m^2)

parameters.CLtable = [
    0.299968350000000 0.343608540000000 0.371341970000000 0.343649450000000 0.319213950000000;
    0.174426140000000 0.205218640000000 0.215918230000000 0.205841410000000 0.176660330000000;
    -0.00943818920000000 -0.0103667400000000 -0.00749612010000000 0.00871683880000000 0.00280499700000000;
    -0.187362540000000 -0.214696040000000 -0.221827340000000 -0.202221580000000 -0.177797020000000;
    -0.300710750000000 -0.344559650000000 -0.365476970000000 -0.352875820000000 -0.302318360000000
];

parameters.CDtable = [
    0.979610970000000 1.13372810000000 1.21510970000000 1.13892950000000 1.01112690000000;
    1.14289370000000 1.34024850000000 1.41052720000000 1.33679060000000 1.14664540000000;
    1.19495980000000 1.45094210000000 1.53724960000000 1.43150070000000 1.21644420000000;
    1.17123540000000 1.33992400000000 1.43592890000000 1.31202820000000 1.12387350000000;
    1.00943130000000 1.13232540000000 1.25727120000000 1.15570860000000 0.996404960000000
];

% Specify the acceptable terminal range error in meters
parameters.RANGE_ERROR_TOLERANCE = 2000;

% Enable/Disable debug plots and statements
parameters.debug = false;

% Calculate the trim aero parameters
[trimAoA, trimCD, trimCL] = calculateTrimAero(parameters, constants, constants.INITIAL_VELOCITY, constants.INITIAL_ALTITUDE);

fprintf('Trim AoA: %.2f degrees, Trim CL: %.4f, Trim CD: %.4f\n', trimAoA, trimCL, trimCD);

% Non-dimensionalization factors
constants.VELOCITY_SCALE = sqrt(constants.GRAVITY_NOMINAL * constants.EARTH_RADIUS_EQ);
constants.LENGTH_SCALE   = constants.EARTH_RADIUS_EQ;
constants.TIME_SCALE     = sqrt(constants.EARTH_RADIUS_EQ / constants.GRAVITY_NOMINAL);

% Simulation parameters
parameters.dt = 0.001;    % Time step (s)
parameters.tMax = 10;     % Maximum simulation time (s)
parameters.iter = 1;      % Number if iterations
parameters.deltaR = 100;  % Perturbation in altitude used for calculating the partial derivative of density w.r.t. r (m)
parameters.guidanceRate = 0.5; % Rate at which guidance is called (s)
parameters.timeTolerance = 1e-4;

% Initialize state variables
[rho, ~, ~] = atmosphereModel(constants.INITIAL_ALTITUDE);
stateHistory = StateHistory();
t = 0;
r = (constants.INITIAL_ALTITUDE + constants.EARTH_RADIUS_EQ) / constants.LENGTH_SCALE;
v = constants.INITIAL_VELOCITY / constants.VELOCITY_SCALE;
s = calculateGCR(constants.INITIAL_LAT, constants.INITIAL_LON, constants.TARGET_LAT, constants.TARGET_LON, constants.EARTH_RADIUS_EQ, 'deg') / constants.LENGTH_SCALE;
state.flightPathAngle        = constants.INITIAL_FPA;
state.rangeToGo              = s;
state.altitudeUnnormalized   = constants.INITIAL_ALTITUDE;
state.velocityUnnormalized   = constants.INITIAL_VELOCITY;
state.rangeToGoUnnormalized  = normalizeState(s,'length',false,constants);
state.normGeocentricDistance = r;
state.velocity               = v;
state.latitude               = deg2rad(constants.INITIAL_LAT);
state.longitude              = deg2rad(constants.INITIAL_LON);
state.heading                = 0;
state.posECEF                = geodeticToECEF(state.latitude, state.longitude, state.altitudeUnnormalized, constants);
state.bankAngle              = constants.INITIAL_BANK_ANGLE;
state.lift                   = 0.5 * rho * trimCL * state.velocityUnnormalized^2 / constants.VEHICLE_WEIGHT;
state.drag                   = 0.5 * rho * trimCD * state.velocityUnnormalized^2 / constants.VEHICLE_WEIGHT;
state.loadFactor             = sqrt(state.lift^2 + state.drag^2);
state.heatingRate            = constants.K_Q * sqrt(rho) * state.velocity^3.15;

% Calculate final energy
finalEnergy = -1/2 * (constants.TARGET_VELOCITY / constants.VELOCITY_SCALE)^2 + ((constants.TARGET_ALTITUDE + constants.EARTH_RADIUS_EQ) / constants.LENGTH_SCALE)^-1;

% Compute current energy
currentEnergy = -1/2 * state.velocity^2 + 1 / state.normGeocentricDistance;

energyStep = (finalEnergy - currentEnergy) / 20000;
lastGuidanceCall = -0.5;
timeHistory = t;
flagFirstPass = true;
energy                   = currentEnergy:energyStep:finalEnergy;    

tic
%% Simulation loop
for i = 1:length(energy)-1
    %t < parameters.tMax && state.altitude > 0
    state.velocity = sqrt(2 * (1 / state.normGeocentricDistance - energy(i)));

    if currentEnergy >= finalEnergy
        break;
    end

    % Predictor-corrector guidance
    if (timeHistory - lastGuidanceCall >= parameters.guidanceRate) || flagFirstPass
        [commandedBankAngle, predictedTrajectory] = predictorCorrectorGuidance(parameters,...
                                                                               constants,...
                                                                               state,...
                                                                               currentEnergy,...
                                                                               finalEnergy);
        state.commandedBankAngle = commandedBankAngle;
        stateHistory = stateHistory.appendState(timeHistory, state);
        lastGuidanceCall = timeHistory;
        sprintf("Guidance being called at time: %.2f", t)
        flagFirstPass = false;
    end

    % Update state variables
    [state, derivatives] = updateStateEnergy(state, commandedBankAngle, constants, parameters, energyStep);

    % Increment time
    t = energyStep / (state.drag * state.velocity);
    timeHistory = t * constants.TIME_SCALE;
    currentEnergy = currentEnergy + energyStep;
    parameters.iter = parameters.iter + 1;
end
toc

disp('fin')

stateHistory.plotState()

% Plot results
% figure;
% plot(t * constants.TIME_SCALE, r * constants.LENGTH_SCALE);
% xlabel('Time (s)');
% ylabel('Altitude (m)');
%
% figure;
% plot(range * constants.LENGTH_SCALE, r * constants.LENGTH_SCALE);
% xlabel('Range (m)');
% ylabel('Altitude (m)');
