clear, close all

% Constants
constants.GRAVITY_NOMINAL       = 9.81;           % Gravitational acceleration (m/s^2)
constants.EARTH_RADIUS_EQ       = 6378137;        % Earth's radius (m)
constants.INITIAL_ALTITUDE      = 122000;         % Initial altitude (m)
constants.INITIAL_VELOCITY      = 10500;          % Initial velocity (m/s)
constants.INITIAL_FPA           = -5.2 * pi/180;  % Initial flight path angle (rad)
constants.INITIAL_LAT           = -4.7;           % EI latitude (deg)
constants.INITIAL_LON           = -112;           % EI longitude (deg)
constants.TARGET_ALTITUDE       = 31000;          % Target final altitude (m)
constants.TARGET_VELOCITY       = 690;            % Target final velocity (m/s)
constants.TARGET_LAT            = 40;             % Final latitude (deg)
constants.TARGET_LON            = -112;           % Final longitude (deg)
constants.VEHICLE_MASS          = 75.7;           % Vehicle mass (kg)
constants.INITIAL_BANK_ANGLE    = deg2rad(0);
constants.FINAL_BANK_ANGLE      = deg2rad(70);
constants.REFERENCE_LENGTH      = 0.69; % reference length for moment calculation (m)
constants.FRONTAL_AREA          = 1;%0.3849; % frontal area of ADEPT craft (m^2)
constants.VEHICLE_WEIGHT        = constants.VEHICLE_MASS * constants.GRAVITY_NOMINAL;
constants.FLATTENING_FACTOR     = 1 / 298.257223563;
constants.EARTH_RADIUS_POLE     = constants.EARTH_RADIUS_EQ * (1 - constants.FLATTENING_FACTOR);
constants.MAX_LOAD_FACTOR       = 3;             % Max load factor (g)
constants.K_Q                   = 9.4369e-5 * (sqrt(constants.GRAVITY_NOMINAL * constants.EARTH_RADIUS_EQ))^3.15;
constants.BANK_ANGLE_RATE_LIMIT = deg2rad(10); % 10 degree bank angle rate limit
constants.MAX_HEATING_RATE      = 2500000; % Maximum heating rate (w/m^2)
constants.C1                    = 5.21e-3 * .5;
constants.C0                    = 8.71e-5 * .5;
constants.EARTH_ROTATION_RATE   = 7.2921151e-5;

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

parameters.CStable =   [
    0.281708410000000	0.163369390000000	0.0	-0.161948400000000	-0.290602000000000;
    0.337056500000000	0.198105980000000	0.0	-0.205414010000000	-0.342958210000000;
    0.352897360000000	0.212576720000000	0.0	-0.215311590000000	-0.388649110000000;
    0.360878870000000	0.209320120000000	0.0	-0.210080210000000	-0.337261610000000;
    0.283326290000000	0.164757430000000	0.0	-0.167250070000000	-0.285702430000000
];

% Specify the acceptable terminal range error in meters
parameters.RANGE_ERROR_TOLERANCE = 0000;

% Enable/Disable debug plots and statements
parameters.debug = false;

% Calculate the trim aero parameters
[trimAoA, trimCD, trimCL] = calculateTrimAero(parameters, constants, constants.INITIAL_VELOCITY, constants.INITIAL_ALTITUDE);

constants.INITIAL_AOA = deg2rad(0);
constants.INITIAL_SIDESLIP = deg2rad(0);
constants.FINAL_ALPHA = deg2rad(20);

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
parameters.guidanceRate = 1; % Rate at which guidance is called (s)
parameters.timeTolerance = 1e-4;
lastGuidanceCall = -0.5;

% Initialize state variables
[rho, ~, ~] = atmosphereModel(constants.INITIAL_ALTITUDE);
[trimCL, trimCD, trimCS] = calculateAeroCoefficients(parameters, constants.INITIAL_AOA, constants.INITIAL_SIDESLIP);
stateHistory = StateHistory();
t = 0;
r = (constants.INITIAL_ALTITUDE + constants.EARTH_RADIUS_EQ) / constants.LENGTH_SCALE;
v = constants.INITIAL_VELOCITY / constants.VELOCITY_SCALE;
s = calculateGCR(constants.INITIAL_LAT, constants.INITIAL_LON, constants.TARGET_LAT, constants.TARGET_LON, constants.EARTH_RADIUS_EQ, 'deg') / constants.LENGTH_SCALE;
[~, terminalAzimuth] = distance('gc',constants.INITIAL_LAT, constants.INITIAL_LON, constants.TARGET_LAT, constants.TARGET_LON);
state.flightPathAngle        = constants.INITIAL_FPA;
state.rangeToGo              = s;
state.altitudeUnnormalized   = constants.INITIAL_ALTITUDE;
state.velocityUnnormalized   = constants.INITIAL_VELOCITY;
state.rangeToGoUnnormalized  = normalizeState(s,'length',false,constants);
state.normGeocentricDistance = r;
state.velocity               = v;
state.latitude               = deg2rad(constants.INITIAL_LAT);
state.longitude              = deg2rad(constants.INITIAL_LON);
state.heading                = deg2rad(terminalAzimuth);
state.posECEF                = lla2ecef([state.latitude, state.longitude, state.altitudeUnnormalized], 'WGS84');
state.bankAngle              = constants.INITIAL_BANK_ANGLE;
state.angleOfAttack          = constants.INITIAL_AOA;
state.sideslipAngle          = constants.INITIAL_SIDESLIP;
state.lift                   = 0.5 * rho * trimCL * state.velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;
state.drag                   = 0.5 * rho * trimCD * state.velocityUnnormalized^2 * constants.FRONTAL_AREA / constants.VEHICLE_WEIGHT;
state.sideForce              = 0;
state.loadFactor             = sqrt(state.lift^2 + state.drag^2);
state.heatingRate            = constants.K_Q * sqrt(rho) * state.velocity^3.15;
state.headingError           = wrapToPi(terminalAzimuth - state.heading);
% [~, state.crossrange]        = calculateDownrangeCrossrange(state,constants);
% state.crossrange             = normalizeState(state.crossrange, 'length', true, constants);
state.crossrange             = calculateCrossrange(state, constants);
state.bankSign               = -sign(state.crossrange + eps);
state.deadbandWidth          = calculateDeadbandWidth(state.velocity, constants);

% Calculate final energy
finalEnergy = -1/2 * (constants.TARGET_VELOCITY / constants.VELOCITY_SCALE)^2 + ((constants.TARGET_ALTITUDE + constants.EARTH_RADIUS_EQ) / constants.LENGTH_SCALE)^-1;

flagInsideDeadband = abs(state.crossrange) < calculateDeadbandWidth(state.velocity, constants);

tic
%% Simulation loop
while state.velocityUnnormalized > constants.TARGET_VELOCITY
    % Compute current energy
    currentEnergy = -1/2 * state.velocity^2 + 1 / state.normGeocentricDistance;

    if currentEnergy >= finalEnergy
        break;
    end

    % Predictor-corrector guidance
    if t == 0 || (t - lastGuidanceCall >= parameters.guidanceRate)
        % Evaluate the bank angle sign based on the crossrange error and deadband
        if state.crossrange > calculateDeadbandWidth(state.velocity, constants)
            state.bankSign = -1;
        elseif state.crossrange < -calculateDeadbandWidth(state.velocity, constants)
            state.bankSign = 1;
        end
        [commandedAngleOfAttack, commandedSideslipAngle, predictedTrajectory] = predictorCorrectorGuidanceURC(parameters,...
                                                                                constants,...
                                                                                state,...
                                                                                currentEnergy,...
                                                                                finalEnergy,...
                                                                                t);
        state.commandedAngleOfAttack = commandedAngleOfAttack;
        state.commandedSideslipAngle = commandedSideslipAngle;
        heatingRate = state.heatingRate;
        stateHistory = stateHistory.appendState(t, state);
        lastGuidanceCall = t;
        sprintf("Guidance being called at time: %.2f\nCrossrange: %.6f\nDeadband Width: %.6f\nBank Sign: %.2f\nBank Angle Command: %.2f", t, state.crossrange, calculateDeadbandWidth(state.velocity, constants), state.bankSign, commandedBankAngle)
    end

    % Update state variables
    [state, derivatives] = updateStateFullURC(state, commandedAngleOfAttack, commandedSideslipAngle, constants, parameters);

    % Increment time

    t = t + parameters.dt;
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

