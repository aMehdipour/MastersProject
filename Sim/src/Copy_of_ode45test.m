clear
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
constants.MU               = 3.986e5;

load std_atmos.mat

% Initialization of parameters
V = 8142;
Vscale = 8142;% / constants.VEL_SCALE;
alt = 100e3;% * constants.LENGTH_SCALE;
V_fscale = 590;% / constants.VEL_SCALE;
alt_f = 40e3;% * constants.LENGTH_SCALE;
t_initial = 0;% / (alt) - (Vscale^2) / (2); % Initial energy state
t_final = 150;% / (alt_f) - (V_fscale^2) / (2); % Final energy state
y_initial = [alt + constants.RADIUS_EQ; deg2rad(-114.2); deg2rad(35); V; deg2rad(-5.5); 0]; % Initial state vector
sigma0 = 0; % Initial bank angle
sigmaf = deg2rad(70); % Final bank angle
% Time span
tspan = [0, t_final]; % start and end time


% Numerical integration using ode45
% [t, y] = ode45(@(t, y) entryDynamics(t, y, sigma0, sigmaf, t_initial, t_final, constants), [t_initial, t_final], y_initial);

% Solve ODEs
[t, y] = ode45(@(t, y) entry_equations(t, y, constants.MU, constants.OMEGA_EARTH, constants, std), tspan, y_initial);


% Plotting the results
for i = 1:length(y_initial)
    if i == 1
        y(:,i) = y(:,i) - constants.RADIUS_EQ;
    end
    figure
    plot(t, y(:,i))
xlabel('time');
ylabel('State variables');
title('Predictor Step Integration Results');
end

% Define the bank-angle profile function
function sigma = bankAngleProfile(t, sigma0, sigmaf, t_initial, t_final)
    sigma = sigma0 + (sigmaf - sigma0) * (t - t_initial) / (t_final - t_initial);
end

function dy = entryDynamics(t, y, sigma0, sigmaf, t_initial, t_final, constants)
    % Extract state variables from y
    r = y(1);% * constants.LENGTH_SCALE;
    lon = y(2);
    lat = y(3);
    V = y(4);
    fpa = y(5);
    psi = y(6);

    % V = sqrt(2 * (1 / r - e));

    % Compute bank angle
    bank = 0;%bankAngleProfile(t, sigma0, sigmaf, t_initial, t_final);

    % Compute aerodynamic forces using the dummy model
    [L, D, S] = ComputeAeroDummy(V, constants);
    S = 0;
    L = L / constants.MASS;
    D = D / constants.MASS;

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
    V_dot = -D - (constants.MU * sinfpa) / r^2 + constants.OMEGA_EARTH^2 * r * coslat * (sinfpa * coslat - cosfpa * sinlat * cospsi);
    fpa_dot = 1 / V * (L * cosbank + (V^2 / r - constants.MU / r^2) * cosfpa + 2 * constants.OMEGA_EARTH * V * coslat * sinpsi);
    psi_dot = 1 / V * ((L * sinbank)/ cosfpa + V^2 / r * cosfpa * sinpsi * tan(lat) - 2 * constants.OMEGA_EARTH * V * (tan(fpa) * cospsi * coslat - sinlat)...
              + (constants.OMEGA_EARTH^2 * r ) / cosfpa * sinpsi * sinlat * coslat);


  %   fpa_dot = 1 / V * (L / constants.MASS * cosbank - S / constants.MASS * sinbank - constants.GRAVITY_NOMINAL * cosfpa + (V^2 / r) * cosfpa + 2 * constants.OMEGA_EARTH * coslat * sinpsi...
  %             + constants.OMEGA_EARTH^2 * r * coslat * (cosfpa * coslat + sinfpa * sinlat * cospsi));
  %   psi_dot = 1 / V * (L / constants.MASS * sinbank / cosfpa + S / constants.MASS * cosbank / cosfpa - V^2 / r * cosfpa * sinpsi * tan(lat) + 2 * constants.OMEGA_EARTH * V *...
  %             (tan(fpa) * coslat * cospsi - sinlat) - (constants.OMEGA_EARTH^2 * r) / cosfpa * sinlat * coslat * sinpsi);
  %   V_dot = -D/constants.MASS - constants.GRAVITY_NOMINAL * sinfpa / r^2 - constants.OMEGA_EARTH^2 * r * coslat * (sinfpa * coslat - cosfpa * sinlat * cospsi);
  %   % fpa_dot = (L * cosbank) / V - (constants.GRAVITY_NOMINAL * cosfpa) / (V * r) - (2 * constants.OMEGA_EARTH * V * coslat * sinpsi) / (V * r) - constants.OMEGA_EARTH^2 * r * coslat * (cosfpa * coslat + sinfpa * cospsi * sinlat);
  %   % psi_dot = (L * sinbank) / (V * cosfpa) + (V * cosfpa * sinpsi * tan(lat)) / r - (2 * constants.OMEGA_EARTH * V * (tan(fpa) * cospsi * coslat + sinlat)) / (V * cosfpa) + (constants.OMEGA_EARTH^2 * r * sinpsi * sinlat * coslat) / (V * cosfpa);
  % 
  % 
  % % r_dot = V * sinfpa;
  %   % lon_dot = (V * cosfpa * sinpsi) / (r * coslat);
  %   % lat_dot = (V * cosfpa * cospsi) / r;
  %   % % V_dot = -D - constants.GRAVITY_NOMINAL * sinfpa / r^2 - constants.OMEGA_EARTH^2 * r * coslat * (sinfpa * coslat - cosfpa * sinlat * cospsi);
  %   % fpa_dot = (L * cosbank) / V - (constants.GRAVITY_NOMINAL * cosfpa) / (V * r) - (2 * constants.OMEGA_EARTH * V * coslat * sinpsi) / (V * r) - constants.OMEGA_EARTH^2 * r * coslat * (cosfpa * coslat + sinfpa * cospsi * sinlat);
  %   % psi_dot = (L * sinbank) / (V * cosfpa) + (V * cosfpa * sinpsi * tan(lat)) / r - (2 * constants.OMEGA_EARTH * V * (tan(fpa) * cospsi * coslat + sinlat)) / (V * cosfpa) + (constants.OMEGA_EARTH^2 * r * sinpsi * sinlat * coslat) / (V * cosfpa);

    % Pack the derivatives into a column vector
    dy = [r_dot; lon_dot; lat_dot; V_dot; fpa_dot; psi_dot];
end


function [L, D] = ComputeAeroDummy(V, rho, constants)
    CL0 = 0.5; % Base coefficient of lift
    % CD0 = 0.05; % Base coefficient of drag
    LDRatio = 4; % Lift-to-drag ratio

    L = CL0 * V^2/2 * rho; % Lift proportional to square of velocity
    D = 1.5 * V^2/2 * rho; % Drag based on lift-to-drag ratio

    S = constants.REFERENCE_LENGTH; % Reference area (square meters, example value)
    L = L * S;
    D = D * S;
end

function dy = entry_equations(t, y, mu, Omega, constants, std)
    % Unpack state variables
    r = y(1);
    theta = y(2);
    phi = y(3);
    V = y(4);
    gamma = y(5);
    psi = y(6);

    alt = r - constants.RADIUS_EQ;

    rho = interp1(std.alt,std.rho,alt,'linear','extrap');

    % Compute aerodynamic forces
    [L, D] = ComputeAeroDummy(V, rho, constants);

    sigma = 0;

    % Equations of motion
    dr = V * sin(gamma);
    dtheta = V * cos(gamma) * sin(psi) / (r * cos(phi));
    dphi = V * cos(gamma) * cos(psi) / r;
    dV = -D - mu * sin(gamma) / r^2 + Omega^2 * r * cos(phi) * (sin(gamma) * cos(phi) - cos(gamma) * sin(phi) * cos(psi));
    dgamma = 0;%1 / V * (L * cos(sigma) + (V^2 / r - mu / r^2) * cos(gamma) + 2 * Omega * V * cos(phi) * sin(psi) + Omega^2 * r * cos(phi) * (cos(gamma) * cos(phi) + sin(gamma) * cos(psi) * sin(phi)));
    dpsi = 1 / V * (L * sin(sigma) / cos(gamma) + V^2 / r * cos(gamma) * sin(psi) * tan(phi) - 2 * Omega * V * (tan(gamma) * cos(psi) * cos(phi) - sin(phi)) + Omega^2 * r * cos(gamma) * sin(psi) * sin(phi) * cos(phi));

     if t > 64
        dummy = 0;
     end

    % Return derivatives
    dy = [dr; dtheta; dphi; dV; dgamma; dpsi];
end