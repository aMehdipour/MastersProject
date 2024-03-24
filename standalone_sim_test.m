clear

% Constants
constants.GRAVITY_NOMINAL  = 9.81;       % Gravitational acceleration (m/s^2) constants.EARTH_RADIUS_EQ = 6378137;    % Earth's radius (m)
constants.INITAL_ALTITUDE  = 100000;     % Initial altitude (m)
constants.INITIAL_VELOCITY = 7000;       % Initial velocity (m/s)
constants.INITIAL_FPA      = -5 * pi/180;  % Initial flight path angle (rad)
constants.TARGET_ALTITUDE  = 50000; % Target final altitude (m)
constants.TARGET_VELOCITY  = 2000;  % Target final velocity (m/s)
constants.VEHICLE_MASS     = 73.8;     % Vehicle mass (kg)
constants.VEHICLE_WEIGHT   = constants.VEHICLE_MASS * constants.GRAVITY_NOMINAL;

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

% Calculate the trim aero parameters
[trimAoA, trimCD, trimCL] = calculateTrimAero(parameters, constants, velocity, altitude);

fprintf('Trim AoA: %.2f degrees, Trim CL: %.4f, Trim CD: %.4f\n', trimAoA, trimCL, trimCD);

% Non-dimensionalization factors
constants.VELOCITY_SCALE = sqrt(constants.GRAVITY_NOMINAL * constants.EARTH_RADIUS_EQ);
constants.LENGTH_SCALE   = constants.EARTH_RADIUS_EQ;
constants.TIME_SCALE     = sqrt(constants.EARTH_RADIUS_EQ / constants.GRAVITY_NOMINAL);

% Simulation parameters
dt = 0.001;    % Time step (s)
tMax = 10;     % Maximum simulation time (s)

% Initialize state variables
t = 0;
r = (constants.INITAL_ALTITUDE + constants.EARTH_RADIUS_EQ) / constants.LENGTH_SCALE;
v = constants.INITIAL_VELOCITY / constants.VELOCITY_SCALE;
gamma = constants.INITIAL_FPA;
range = 3704000;

% Calculate final energy
finalEnergy = -1/2 * (constants.TARGET_VELOCITY / constants.VELOCITY_SCALE)^2 + ((constants.TARGET_ALTITUDE + constants.EARTH_RADIUS_EQ) / constants.LENGTH_SCALE)^-1;

% Initialize guidance parameters
bankAngle = 0;

% Simulation loop
while t < tMax && r > 0
    % Compute current energy
    currentEnergy = -1/2 * v^2 + 1 / r;
    
    % Predictor-corrector guidance
    [bankAngle, ~] = predictorCorrectorGuidance(currentEnergy, finalEnergy, range, trimCL, trimCD);
    
    % Update state variables
    [r, v, gamma, range] = updateState(r, v, gamma, range, bankAngle, dt, trimCL, trimCD);
    
    % Increment time
    t = t + dt;
end

% Plot results
figure;
plot(t * constants.TIME_SCALE, r * constants.LENGTH_SCALE);
xlabel('Time (s)');
ylabel('Altitude (m)');

figure;
plot(range * constants.LENGTH_SCALE, r * constants.LENGTH_SCALE);
xlabel('Range (m)');
ylabel('Altitude (m)');



% Evaluate terminal constraints function
function terminalConstraint = evaluateTerminalConstraints(predictedTrajectory, finalEnergy)
    % Extract final energy from predicted trajectory
    finalEnergy_predicted = predictedTrajectory.energy(end);
    
    % Compute terminal constraint
    terminalConstraint = finalEnergy_predicted - finalEnergy;
end

% Update state variables function
function [altitude, velocity, flightPathAngle, range] = updateState(altitude, velocity, flightPathAngle, range, bankAngle, dt, CL, CD)
    % Compute non-dimensionalized lift and drag forces
    [rho, ~, ~] = atmosphereModel(altitude);
    lift = 0.5 * rho * CL * velocity^2;
    drag = 0.5 * rho * CD * velocity^2;
    
    % Equations of motion
    velocityDot = -drag - sin(flightPathAngle);
    flightPathAngleDot = (lift * cos(bankAngle) - cos(flightPathAngle)) / velocity;
    altitudeDot = velocity * sin(flightPathAngle);
    rangeDot = velocity * cos(flightPathAngle);
    
    % Update state variables
    velocity = velocity + velocityDot * dt;
    flightPathAngle = flightPathAngle + flightPathAngleDot * dt;
    altitude = altitude + altitudeDot * dt;
    range = range + rangeDot * dt;
end
