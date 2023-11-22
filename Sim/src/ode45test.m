%% Constants
c = 271.42; % speed of sound at 100km altitude (m/s)
constants.RADIUS_EQ        = 6378137.0;    % Earth's radius at the equator (semi-major axis)
constants.RADIUS_POLE      = 6356752.3142; % Earth's radius at the poles (semi-minor axis)
constants.MASS             = 73.8;%8.490; % vehicle mass (kg)
constants.INERTIA          = [0.25 0 0; 0 0.1722 0; 0 0 0.1719]; % vehicle moment of inertia tensor (kg-m^2)
constants.FLAP_AREA        = 0.026;
constants.FRONTAL_AREA     = 0.3849; % frontal area of ADEPT craft (m^2)
constants.REFERENCE_LENGTH = 0.69; % reference length for moment calculation (m)
constants.FLAP_MOMENT_ARM  = 0.84;
constants.FLAP_LIMIT       = 20; % Flap deflection limit (degrees)
constants.OMEGA_EARTH      = 7.292115e-5; % Earth's rotation rate (rad/s)
constants.GRAVITY_NOMINAL  = 9.81; % Nominal gravity at sea level (m/s^2)
constants.TIME_SCALE       = sqrt(constants.RADIUS_EQ / constants.GRAVITY_NOMINAL);
constants.VEL_SCALE        = sqrt(constants.GRAVITY_NOMINAL * constants.RADIUS_EQ);
constants.LENGTH_SCALE     = 1 / constants.RADIUS_EQ;
constants.ALTITUDE_FINAL   = 6000; % Target altitude for guidance termination
constants.RANGE_FINAL      = 2000; % Allowable miss distance
constants.VELOCITY_FINAL   = 690; % Velocity magnitude at which to exit guidance (m/s). Taken from 
constants.G_LIMIT          = 15; % Structural G limit in gees
constants.HEAT_RATE_LIMIT  = 250; % Heat rate limit in W/cm^2

% Initialization of parameters
V = 8142;
Vscale = 8142;% / constants.VEL_SCALE;
alt = 100e3;% * constants.LENGTH_SCALE;
V_fscale = 590;% / constants.VEL_SCALE;
alt_f = 40e3;% * constants.LENGTH_SCALE;
e_initial = 1 / (alt) - (Vscale^2) / (2); % Initial energy state
e_final = 1 / (alt_f) - (V_fscale^2) / (2); % Final energy state
y_initial = [alt; deg2rad(-114.2); deg2rad(35); deg2rad(-5.5); 0]; % Initial state vector
sigma0 = 0; % Initial bank angle
sigmaf = -pi; % Final bank angle


% Numerical integration using ode45
[e, y] = ode45(@(e, y) entryDynamics(e, y, sigma0, sigmaf, e_initial, e_final, constants), [e_initial, e_final], y_initial);

% Plotting the results
figure
plot(e, y);
legend('State 1', 'State 2', 'State 3');
xlabel('Energy-like state (e)');
ylabel('State variables');
title('Predictor Step Integration Results');

% Define the bank-angle profile function
function sigma = bankAngleProfile(e, sigma0, sigmaf, e_initial, e_final)
    sigma = sigma0 + (sigmaf - sigma0) * (e - e_initial) / (e_final - e_initial);
end

function dy = entryDynamics(e, y, sigma0, sigmaf, e_initial, e_final, constants)
    % Extract state variables from y
    r = y(1);% * constants.LENGTH_SCALE;
    lon = y(2);
    lat = y(3);
    % V = y(4);
    fpa = y(4);
    psi = y(5);

    V = sqrt(2 * (1 / r - e));

    % Compute bank angle
    bank = bankAngleProfile(e, sigma0, sigmaf, e_initial, e_final);

    % Compute aerodynamic forces using the dummy model
    [L, D, S] = ComputeAeroDummy(V, constants);
    S = 0;

    % Pre-compute trigonometric values
    sinfpa = sin(fpa);
    cosfpa = cos(fpa);
    sinpsi = sin(psi);
    cospsi = cos(psi);
    sinlat = sin(lat);
    coslat = cos(lat);
    cosbank = cos(bank);
    sinbank = sin(bank);

    % EOMs
    r_dot = V * sinfpa;
    lon_dot = (V * cosfpa * sinpsi) / (r * coslat);
    lat_dot = (V * cosfpa * cospsi) / r;
    fpa_dot = 1 / V * (L / constants.MASS * cosbank - S / constants.MASS * sinbank - constants.GRAVITY_NOMINAL * cosfpa + (V^2 / r) * cosfpa + 2 * constants.OMEGA_EARTH * coslat * sinpsi...
              + constants.OMEGA_EARTH^2 * r * coslat * (cosfpa * coslat + sinfpa * sinlat * cospsi));
    psi_dot = 1 / V * (L / constants.MASS * sinbank / cosfpa + S / constants.MASS * cosbank / cosfpa - V^2 / r * cosfpa * sinpsi * tan(lat) + 2 * constants.OMEGA_EARTH * V *...
              (tan(fpa) * coslat * cospsi - sinlat) - (constants.OMEGA_EARTH^2 * r) / cosfpa * sinlat * coslat * sinpsi);
    % % V_dot = -D - constants.GRAVITY_NOMINAL * sinfpa / r^2 - constants.OMEGA_EARTH^2 * r * coslat * (sinfpa * coslat - cosfpa * sinlat * cospsi);
    % fpa_dot = (L * cosbank) / V - (constants.GRAVITY_NOMINAL * cosfpa) / (V * r) - (2 * constants.OMEGA_EARTH * V * coslat * sinpsi) / (V * r) - constants.OMEGA_EARTH^2 * r * coslat * (cosfpa * coslat + sinfpa * cospsi * sinlat);
    % psi_dot = (L * sinbank) / (V * cosfpa) + (V * cosfpa * sinpsi * tan(lat)) / r - (2 * constants.OMEGA_EARTH * V * (tan(fpa) * cospsi * coslat + sinlat)) / (V * cosfpa) + (constants.OMEGA_EARTH^2 * r * sinpsi * sinlat * coslat) / (V * cosfpa);


  % r_dot = V * sinfpa;
    % lon_dot = (V * cosfpa * sinpsi) / (r * coslat);
    % lat_dot = (V * cosfpa * cospsi) / r;
    % % V_dot = -D - constants.GRAVITY_NOMINAL * sinfpa / r^2 - constants.OMEGA_EARTH^2 * r * coslat * (sinfpa * coslat - cosfpa * sinlat * cospsi);
    % fpa_dot = (L * cosbank) / V - (constants.GRAVITY_NOMINAL * cosfpa) / (V * r) - (2 * constants.OMEGA_EARTH * V * coslat * sinpsi) / (V * r) - constants.OMEGA_EARTH^2 * r * coslat * (cosfpa * coslat + sinfpa * cospsi * sinlat);
    % psi_dot = (L * sinbank) / (V * cosfpa) + (V * cosfpa * sinpsi * tan(lat)) / r - (2 * constants.OMEGA_EARTH * V * (tan(fpa) * cospsi * coslat + sinlat)) / (V * cosfpa) + (constants.OMEGA_EARTH^2 * r * sinpsi * sinlat * coslat) / (V * cosfpa);

    % Pack the derivatives into a column vector
    dy = [r_dot; lon_dot; lat_dot; fpa_dot; psi_dot];
end


function [L, D, S] = ComputeAeroDummy(V, constants)
    CL0 = 0.5; % Base coefficient of lift
    % CD0 = 0.05; % Base coefficient of drag
    LDRatio = 4; % Lift-to-drag ratio

    L = CL0 * V^2; % Lift proportional to square of velocity
    D = L / LDRatio; % Drag based on lift-to-drag ratio

    S = 1; % Reference area (square meters, example value)
    L = L * S;
    D = D * S;
end
