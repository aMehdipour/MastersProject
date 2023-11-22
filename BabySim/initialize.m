function [constants, state, sim] = intitialize()

  %% Inputs and constants
  constants.c                = 271.42; % speed of sound at 100km altitude (m/s)
  constants.INITIAL_MACH     = 30; % Mach number (dless) and speed of sound (m/s)
  constants.INITIAL_ALT      = 100e3;
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
  constants.TIME_SCALE       = 1 / sqrt(constants.RADIUS_EQ / constants.GRAVITY_NOMINAL);
  constants.VEL_SCALE        = 1 / sqrt(constants.GRAVITY_NOMINAL * constants.RADIUS_EQ);
  constants.LENGTH_SCALE     = 1 / constants.RADIUS_EQ;
  constants.ALTITUDE_FINAL   = 40e3; % Target altitude for guidance termination
  constants.RANGE_FINAL      = 5000; % Allowable miss distance
  constants.VELOCITY_FINAL   = 690; % Velocity magnitude at which to exit guidance (m/s). Taken from 
  constants.G_LIMIT          = 15; % Structural G limit in gees
  constants.HEAT_RATE_LIMIT  = 250; % Heat rate limit in W/cm^2
  constants.LAT_TARGET       = deg2rad(37.3352); % SJSU latitude 37.3352 degrees
  constants.LON_TARGET       = deg2rad(121.8811); % SJSU longitude 121.8811 degrees

  %% Sim Settings
  sim.dt = 1/500;
  sim.tend = 200;
  sim.tvec = 0:sim.dt:sim.tend;
  sim.npoints = length(sim.tvec); % Number of timesteps

  %% Set up state variables
  % Dynamics
  state.V               = zeros(sim.npoints,1);
  state.Vdot            = state.V;% vehicle velocity (m/s)
  state.flightPathAngle = zeros(sim.npoints,1);
  state.gammaDot        = state.flightPathAngle; % flight path angle (rad)
  state.heading         = zeros(sim.npoints,1);
  state.headingDot      = state.heading; % headinging angle (rad)
  state.chi             = zeros(sim.npoints,1);
  state.chidot          = state.chi; % track angle (rad)
  state.p               = zeros(sim.npoints,1);
  state.pDot            = state.p; % angular velocity about x axis (rad/s)
  state.q               = zeros(sim.npoints,1);
  state.qDot            = state.q; % angular velocity about y axis (rad/s)
  state.r               = zeros(sim.npoints,1);
  state.rDot            = state.r; % angular velocity about z axis (rad/s)
  state.pDotFilt        = state.r;
  state.qDotFilt        = state.r;
  state.rDotFilt        = state.r;
  state.gammaDotFilt    = state.r;
  state.headingDotFilt  = state.r;
  state.rotBodyFromECEF = zeros(3,3,sim.npoints);
  state.rotNedFromBody  = zeros(3,3,sim.npoints);
  state.rotNedFromECEF  = zeros(3,3,sim.npoints);

  % Kinematics
  state.z      = zeros(sim.npoints,1);
  state.zdot   = state.z; % orbital radius (m)
  state.lat    = zeros(sim.npoints,1);
  state.latdot = state.lat; % latitude (rad)
  state.lon    = zeros(sim.npoints,1);
  state.londot = state.lon; % longitude (rad)
  state.bank   = zeros(sim.npoints,1);
  state.bdot   = state.bank; % bank angle (rad)
  state.aoa    = zeros(sim.npoints,1);
  state.adot   = state.aoa; % angle of attack (rad)
  state.ss     = zeros(sim.npoints,1);
  state.ssdot  = state.ss; % sideslip angle (rad)
  state.pitch  = zeros(sim.npoints,1);
  state.roll   = zeros(sim.npoints,1);
  state.alt    = state.z;

  %% Initialize state variables
  state.V(1)       = constants.INITIAL_MACH * constants.c;
  state.heading(1) = 0;
  state.p(1)       = 0; 
  state.q(1)       = 0;
  state.r(1)       = 0;
  state.alt(1)     = constants.INITIAL_ALT;
  state.z(1)       = state.alt(1) + constants.RADIUS_EQ;

  % Set up energy terms
  sim.e_initial = 1 / (constants.INITIAL_ALT) - ((state.V * constants.VEL_SCALE)^2) / (2); % Initial energy state
  sim.e_final = 1 / (constants.INITIAL_ALT * constants.LENGTH_SCALE) - ((constants.VELOCITY_FINAL * constants.VEL_SCALE)^2) / (2); % Final energy state
