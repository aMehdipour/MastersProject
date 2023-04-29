%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADEPTSim
% This script simulates the atmospheric entry of an ADEPT-class entry
% vehicle utilizing an aerodynamic flap control system.
%
% Created By: Joshua Stokes
% Last Modified: 8/14/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; %close all
addpath('/home/amehdipour/repos/MastersProject/Sim/Data');
addpath('/home/amehdipour/repos/MastersProject/Sim/QCAT/qcat');
addpath('/home/amehdipour/repos/MastersProject/Sim/src');
% addpath Outputs
load CFD; load std_atmos; load att_profile
set(0,'DefaultFigureWindowStyle','docked')

%% Constants
c = 271.42; % speed of sound at 100km altitude (m/s)
constants.RADIUS_EQ        = 6378137.0;    % Earth's radius at the equator (semi-major axis)
constants.RADIUS_POLE      = 6356752.3142; % Earth's radius at the poles (semi-minor axis)
constants.MASS             = 8.490; % vehicle mass (kg)
constants.INERTIA          = [0.25 0 0; 0 0.1722 0; 0 0 0.1719]; % vehicle moment of inertia tensor (kg-m^2)
constants.FLAP_AREA        = 0.026;
constants.FRONTAL_AREA     = 0.3849; % frontal area of ADEPT craft (m^2)
constants.REFERENCE_LENGTH = 0.69; % reference length for moment calculation (m)
constants.FLAP_MOMENT_ARM  = 0.84;
constants.FLAP_LIMIT       = 20; % Flap deflection limit (degrees)
constants.OMEGA_EARTH      = 7.292115e-5; % Earth's rotation rate (rad/s)
rng(69); % Set the RNG seed for reproducability

%% Sim Settings
dt = 1/500;
tend = 20;
tvec = 0:dt:tend;
npoints = length(tvec); % Number of timesteps
CFD.aeroScalingFactor = 1; % Scaling factor for aero variables to add scaling error

%% Controller
ctrl.ts = 1/500; % sample time (s)
ts_g = 1/500;
ctrl.tsratio = ts_g/ctrl.ts;
ctrl.cnt = 1;
methods = [0, 1]; % 0 to use NDI,1 to use INDI
ctrl.loopselect = 0; % 0 to use guidance, 1 to use angle, 2 to use rate
cs_rate = 200; % maximum deflection rate for control surfaces (deg/s)
ctrl.maxdeflect = ctrl.ts * cs_rate; %maximum deflection over one ctrl sample (deg)
ctrl.ratelim = 10 * pi/180;
ctrl.anglim = 18 * pi/180;

% Rotational rate controller
ctrl.kp_r = 18;

% Angular displacement controller
ctrl.kp_a = 6;

% Guidance Control
ctrl.kp_g = 1.92;
ctrl.kd_g = 0;
ctrl.ki_g = 0.4;

% Transfer functions (analysis purposes only)
DI = tf([0 1],[1 0]); % dynamic inversion modeled as simple integrator
OL_rate = ctrl.kp_r * DI;
CL_rate = feedback(OL_rate,1); % closed rate loop

OL_angle = ctrl.kp_a * CL_rate;
CL_angle = feedback(OL_angle,1);

% pidTuner(CL_rate)
% pidTuner(DI * CL_angle)
CNT_guidance = tf([ctrl.kd_g ctrl.kp_g ctrl.ki_g],[0 1 0]);

% Active set inputs
ctrl.Wv = [1 0 0; 0 1 0; 0 0 1];%eye(3);
ctrl.Wu = eye(8);
ctrl.gam_wls = 1e6;


for method = methods
    flagFirstPass1 = 0;
    flagFirstPass2 = 0;
    if method == 1 % implement filter for INDI

        % Angular accel filter
        [b,a] = butter(1,0.5);
        ctrl.filt_disc = tf(b,a,ctrl.ts);
        a = [a 0]; b = [b 0];
        ctrl.a = a; ctrl.b = b;

        % figure
        % bode(ctrl.filt_disc)
        % title('Angular Rate Filter')

        % Guidance filter
        [b_g,a_g] = butter(1,250/500);
        ctrl.filt_disc = tf(b_g,a_g,ctrl.ts);
        a_g = [a_g 0]; b_g = [b_g 0];
        ctrl.a_g = a_g; ctrl.b_g = b_g;

        % figure
        % opts1=bodeoptions('cstprefs');
        % opts1.PhaseVisible = 'off';
        % opts1.YLim={[-30 0]};
        % Mag=subplot(2,1,1);bodeplot(ctrl.filt_disc,opts1); grid on;
        % title('INDI Filter')
        % set(xlabel(''),'visible','off');
        % opts2=bodeoptions('cstprefs');
        % opts2.MagVisible = 'off';
        % opts2.YLim={[-90 0]};
        % Phase=subplot(2,1,2);bodeplot(ctrl.filt_disc,opts2); grid on; title('');
    else % set filter as pass-through
        a = [1 0 0]; b = [1 0 0];
        ctrl.a = a; ctrl.b = b;
        a_g = [1 0 0]; b_g = [1 0 0];
        ctrl.a_g = a_g; ctrl.b_g = b_g;
    end

    % Initialize filter memory to zero
    % State feedback filters
    p_in = zeros(1,3); q_in = zeros(1,3); r_in = zeros(1,3); %input memory
    g_in = zeros(1,3); h_in = zeros(1,3);
    p_mem = zeros(1,2); q_mem = zeros(1,2); r_mem = zeros(1,2); %output memory
    g_mem = zeros(1,2); h_mem = zeros(1,2);

    % Control signal filters
    ctrl.vmem = zeros(3,2); %input memory
    ctrl.gmem = zeros(2,2);
    ctrl.fmem = zeros(3,2); %output memory
    ctrl.gfmem = zeros(2,2);

    %% Set up state variables
    % Dynamics
    V = zeros(npoints,1); Vdot = V;% vehicle velocity (m/s)
    flightPathAngle = zeros(npoints,1); gammaDot = flightPathAngle; % flight path angle (rad)
    heading = zeros(npoints,1); headingDot = heading; % headinging angle (rad)
    chi = zeros(npoints,1); chidot = chi; % track angle (rad)
    p = zeros(npoints,1); pDot = p; % angular velocity about x axis (rad/s)
    q = zeros(npoints,1); qDot = q; % angular velocity about y axis (rad/s)
    r = zeros(npoints,1); rDot = r; % angular velocity about z axis (rad/s)
    pDotFilt = r; qDotFilt = r; rDotFilt = r;
    gammaDotFilt = r; headingDotFilt = r;
    rotBodyFromECEF = zeros(3,3,npoints);
    rotNedFromBody  = zeros(3,3,npoints);
    rotNedFromECEF  = zeros(3,3,npoints);

    % Kinematics
    z     = zeros(npoints,1); zdot = z; % orbital radius (m)
    lat   = zeros(npoints,1); latdot = lat; % latitude (rad)
    lon   = zeros(npoints,1); londot = lon; % longitude (rad)
    bank  = zeros(npoints,1); bdot = bank; % bank angle (rad)
    aoa   = zeros(npoints,1); adot = aoa; % angle of attack (rad)
    ss    = zeros(npoints,1); ssdot = ss; % sideslip angle (rad)
    pitch = zeros(npoints,1);
    roll  = zeros(npoints,1);

    % Telemetry
    vlin = zeros(3,npoints-1);
    flin = zeros(3,npoints-1);
    if ~exist('cmd_a'); cmd_a = zeros(3,npoints); end
    cmd_r = zeros(3,npoints-1);
    cmd_g = zeros(2,npoints-1);

    %% Initial Conditions
    M = 30; % Mach number (dless) and speed of sound (m/s)
    V(1) = M * c;
    heading(1) = 0;
    p(1) = 0; q(1) = 0; r(1) = 0;
    alt = 100000; %118000;
    alt_init = alt; % save initial altitude for filename
    z(1) = alt + constants.RADIUS_EQ;

    % SJSU latitude 37.3352 degrees
    % SJSU longitude 121.8811 degrees
    lat(1) = deg2rad(37.3352); lon(1) = deg2rad(121.8811);
    bank(1) = 0; aoa(1) = 0; ss(1) = 0;

    flapDeflection = zeros(8,npoints); % flap deflection angles (deg)
    fc =zeros(8,1);

    % controller states
    ctrl.errorLastGammaHeading = [0; 0];
    ctrl.errorIntegGammaHeading = [0; 0];
    ctrl.fclast = fc;

    %% Case Specific Overrides
    % This section overrides the vehicle initial conditions to result in trim
    % conditions at the start of the simulation

    if 0 %M30 A100
        fc = [-0.7761   -0.9669   -1.4826    1.1879    1.5912    2.6582    0.3138   -0.0417]';
    elseif 0 %M30 A100 FPA-5.5
        aoa(1) = 2.395 * pi/180; ss(1) = -0.400 * pi/180;
        fc = [1.6441    0.1347   -1.6763   -1.3564   -1.2157    0.5829    1.3549    1.8832]';
    end

    % IMU gyro noise specifications
    gammas = -5.5;%-40:20:40;
    gyroNoiseDensity = gammas;%[0, 1e-5, 1e-4, 1e-3];%linspace(0, 1e-3, 10);
    dataSave.gyroNoiseDensity = gyroNoiseDensity;

    %% Execute Computations
    tic
    for j =  1:length(gyroNoiseDensity)
        % V = zeros(npoints,1); Vdot = V;% vehicle velocity (m/s)
        % flightPathAngle = zeros(npoints,1); gammaDot = flightPathAngle; % flight path angle (rad)
        % heading = zeros(npoints,1); headingDot = heading; % headinging angle (rad)
        % chi = zeros(npoints,1); chidot = chi; % track angle (rad)
        % p = zeros(npoints,1); pDot = p; % angular velocity about x axis (rad/s)
        % q = zeros(npoints,1); qDot = q; % angular velocity about y axis (rad/s)
        % r = zeros(npoints,1); rDot = r; % angular velocity about z axis (rad/s)
        % pDotFilt = r; qDotFilt = r; rDotFilt = r;
        % gammaDotFilt = r; headingDotFilt = r;
        % rotBodyFromECEF = zeros(3,3,npoints);
        % rotNedFromBody  = zeros(3,3,npoints);
        % rotNedFromECEF  = zeros(3,3,npoints);

        % % Kinematics
        % z     = zeros(npoints,1); zdot = z; % orbital radius (m)
        % lat   = zeros(npoints,1); latdot = lat; % latitude (rad)
        % lon   = zeros(npoints,1); londot = lon; % longitude (rad)
        % bank  = zeros(npoints,1); bdot = bank; % bank angle (rad)
        % aoa   = zeros(npoints,1); adot = aoa; % angle of attack (rad)
        % ss    = zeros(npoints,1); ssdot = ss; % sideslip angle (rad)
        % pitch = zeros(npoints,1);
        % roll  = zeros(npoints,1);

        % % Telemetry
        % vlin = zeros(3,npoints-1);
        % flin = zeros(3,npoints-1);
        % if ~exist('cmd_a'); cmd_a = zeros(3,npoints); end
        % cmd_r = zeros(3,npoints-1); % [p, q, r]
        % cmd_g = zeros(2,npoints-1); % [flightPathAngle; heading]

        % %% Initial Conditions
        % M = 30; % Mach number (dless) and speed of sound (m/s)
        % V(1) = M * c;
        % heading(1) = 0;
        % p(1) = 0; q(1) = 0; r(1) = 0;
        % alt = 100000; %118000;
        % alt_init = alt; % save initial altitude for filename
        % z(1) = alt + constants.RADIUS_EQ;
        % lat(1) = deg2rad(37.3352); lon(1) = deg2rad(121.8811);
        % bank(1) = 0; aoa(1) = 0; ss(1) = 0;

        % flapDeflection = zeros(8,npoints); % flap deflection angles (deg)
        % fc = zeros(8,1);

        % controller states
        ctrl.errorLastGammaHeading = [0; 0];
        ctrl.errorIntegGammaHeading = [0; 0];
        ctrl.fclast = fc;
        flightPathAngle(1) = gammas(j) * pi/180;
        pitch(1) = flightPathAngle(1) + aoa(1);
        roll(1)  = bank(1);
        chi(1)   = calculateTrackAngle(V(1), lat(1), lon(1), heading(1), ss(1), constants, dt);

        % Rotation matrices
        rotNedFromBody(:,:,1)  = angle2dcm(heading(1),pitch(1),roll(1),'ZYX');
        rotNedFromECEF(:,:,1)  = calcRotNedFromEcef(lat(1),lon(1));
        rotBodyFromECEF(:,:,1) = rotNedFromBody(:,:,1)' * rotNedFromECEF(:,:,1);

        for i = 1:npoints-1 % npoints is the number of time steps

            % Interpolate density and gravitational force
            alt = z(i) - constants.RADIUS_EQ; % TODO: Currently referencing Earth
            rho = interp1(std.alt,std.rho,alt,'linear','extrap'); % atmospheric density (kg/m^3)
            gz = constants.MASS * interp1(std.alt,std.g,alt,'linear','extrap'); % gravitational force (N)

            % Execute Controller
            if mod(i * dt,ctrl.ts) == 0 && i >1

                % Compute aerodynamic force/moment inputs to controller
                [v_aero,B,bias,f_aero,Ffit] = ControllerAero(V(i),rho,aoa(i),ss(i),bank(i),flapDeflection(:,i-1),CFD,constants);

                % Save state variables into input structure
                st = struct('flightPathAngle',flightPathAngle(i),'V',V(i),'pitch',pitch(i),'roll',roll(i),'heading',heading(i),'p',p(i),'q',q(i),'r',r(i),'bank',bank(i),'rho',rho,'flapDeflection',flapDeflection(:,i-1),...
                    'gammaDot',gammaDotFilt(i-1),'z',z(i-1),'headingDot',headingDotFilt(i-1),'pDot',pDotFilt(i-1),'qDot',qDotFilt(i-1),'rDot',rDotFilt(i-1),...
                    'aoa',aoa(i),'ssl',ss(i),'gz',gz,'latitude',lat(i),'longitude',lon(i),'rotBodyFromECEF',rotBodyFromECEF(:,:,i),'rotNedFromBody',rotNedFromBody(:,:,i),...
                    'rotNedFromECEF',rotNedFromECEF(:,:,i));

                % Commands (Input selected based on ctrl.loopselect)
                % NOTE: The original implementation only set the flightPathAngle
                % command to the inital angle for all time.
                if i * dt >= 2
                    cmd_slope = deg2rad(0.02);
                    ctrl.cmd_g = [flightPathAngle(1)-cmd_slope * i * dt + 2 * cmd_slope; 0];
                    if flagFirstPass1 == 0
                        sprintf('The flight path angle command is %.3f', rad2deg(ctrl.cmd_g(1)))
                        flagFirstPass1 = flagFirstPass1 + 1;
                    end
                else
                    ctrl.cmd_g = [flightPathAngle(1) ; 0];
                    if flagFirstPass2 == 0
                        sprintf('The flight path angle command is %.3f', rad2deg(ctrl.cmd_g(1)))
                        flagFirstPass2 = flagFirstPass2 + 1;
                    end
                end

                %         cmd_g(1,1) = -5.5 * pi/180;
                %         ctrl.cmd_g = [cmd_g(1,i-1);0] + [0.75 * V(i)/z(i) * cos(flightPathAngle(i));0] * dt;
                %         cmd_slope = 0.02 * pi/180; %rad/sec
                %         if i * dt < 10
                %             ctrl.cmd_g = [flightPathAngle(1); heading(1)] - cmd_slope * [i * dt; 0];
                %         else
                %             ctrl.cmd_g = [cmd_g(1,i-1); heading(1)] + cmd_slope * [0; i * dt-10];
                %         end
                %         if i == 2
                %             ctrl.cmd_a = [0; aoa(2); ss(2)];
                %         end

                % TODO: Why are these if statements giving commands at specific times?
                if ctrl.loopselect == 1 % If controller is in Angle loop
                    if i * dt > 1 && i * dt < 3 % If the time is between 2 and 3 seconds
                        ctrl.cmd_a = [0; 1; 0] * pi/180; % cmd_a = [bank; aoa; ss]
                    elseif i * dt > 1
                        ctrl.cmd_a = [0; 1; 1] * pi/180;
                    else
                        ctrl.cmd_a = [0; 0; 0] * pi/180;
                    end
                    %               ctrl.cmd_a = cmd_a(:,i);
                end

                % if i * dt > 0.5 && i * dt < 1.5
                %     ctrl.cmd_r = [0; 1; 0] * pi/180;
                % elseif i * dt > 0.5
                %     ctrl.cmd_r = [0; 1; 1] * pi/180;
                % else
                %     ctrl.cmd_r = [0; 0; 0];
                % end

                % Execute control law
                st_sensor = st;
                % st_sensor.q = st_sensor.q + randn * sqrt(gyroNoiseDensity(j));
                % st_sensor.p = st_sensor.p + randn * sqrt(gyroNoiseDensity(j));
                % st_sensor.r = st_sensor.r + randn * sqrt(gyroNoiseDensity(j));
                [fc,ctrl] = RunController(st_sensor,ctrl,v_aero,f_aero,Ffit,B,bias,method,constants);

                %store telemetry
                vlin(:,i) = ctrl.vlin;
                flin(:,i) = ctrl.flin;
                cmd_a(:,i) = ctrl.cmd_a;
                cmd_r(:,i) = ctrl.cmd_r;
                cmd_g(:,i) = ctrl.cmd_g;
            end % end if mod(i * dt,ctrl.ts) == 0

            % Update flap Positions
            flapDeflection(:,i) = fc;

            % Compute aerodynamic forces and moments
            [L,D,S,Lm,Mm,Nm] = ComputeAero(V(i),rho,aoa(i),ss(i),bank(i),flapDeflection(:,i),CFD,constants);

            % Compute State Derivatives
            Vdot(i) = -(gz/constants.MASS) * sin(flightPathAngle(i)) - (1/constants.MASS) * D;
            gammaDot(i) = (-gz/(constants.MASS * V(i)) + (V(i)/z(i))) * cos(flightPathAngle(i)) + 1/(constants.MASS * V(i)) * (L * cos(bank(i)) - S * sin(bank(i)));
            headingDot(i) = (1/(constants.MASS * V(i) * cos(flightPathAngle(i)))) * (L * sin(bank(i)) + S * cos(bank(i)));

            gyro = cross([p(i); q(i); r(i)],constants.INERTIA * [p(i); q(i); r(i)]); %gyroscopic torque (Nm)
            Omegadot = inv(constants.INERTIA) * ([Lm; Mm; Nm;] - gyro); % rotational dynamics derivatives
            pDot(i) = Omegadot(1);
            qDot(i) = Omegadot(2);
            rDot(i) = Omegadot(3);
            dataSave.pDotSave(i, j, method + 1) = pDot(i);
            dataSave.qDotSave(i, j, method + 1) = qDot(i);
            dataSave.rDotSave(i, j, method + 1) = rDot(i);

            zdot(i) = V(i) * sin(flightPathAngle(i));
            bdot(i) = p(i) * cos(aoa(i)) * sec(ss(i)) + r(i) * sin(aoa(i)) * sec(ss(i));
            adot(i) = -p(i) * cos(aoa(i)) * tan(ss(i)) + q(i) - r(i) * sin(aoa(i)) * tan(ss(i));
            ssdot(i) = p(i) * sin(aoa(i)) - r(i) * cos(aoa(i));

            %Filter derivatives for INDI processing
            % Store filter inputs [in(n) in(n-1) in(n-2)]
            p_in = [pDot(i) p_in(1:2)]; q_in = [qDot(i) q_in(1:2)]; r_in = [rDot(i) r_in(1:2)];
            g_in = [gammaDot(i) g_in(1:2)]; h_in = [headingDot(i) h_in(1:2)];
            % Execute filter function
            pDotFilt(i) = filter_discrete(p_in, p_mem,a,b);
            qDotFilt(i) = filter_discrete(q_in, q_mem,a,b);
            rDotFilt(i) = filter_discrete(r_in, r_mem,a,b);
            gammaDotFilt(i) = filter_discrete(g_in, g_mem,a_g,b_g);
            headingDotFilt(i) = filter_discrete(h_in, h_mem,a_g,b_g);
            dataSave.pDotFilt(i, j, method + 1) = rad2deg(pDotFilt(i));
            dataSave.qDotFilt(i, j, method + 1) = rad2deg(qDotFilt(i));
            dataSave.rDotFilt(i, j, method + 1) = rad2deg(rDotFilt(i));

            % Store filter outputs into "memory" [out(n) out(n-1)]
            p_mem = [pDotFilt(i) p_mem(1)]; q_mem = [qDotFilt(i) q_mem(1)]; r_mem = [rDotFilt(i) r_mem(1)];
            g_mem = [gammaDotFilt(i) g_mem(1)]; h_mem = [headingDotFilt(i) h_mem(1)];

            % Update State Variables
            V(i+1) = V(i) + Vdot(i) * dt;
            flightPathAngle(i+1) = flightPathAngle(i) + gammaDot(i) * dt;
            heading(i+1) = heading(i) + headingDot(i) * dt;
            p(i+1) = p(i) + pDot(i) * dt;
            q(i+1) = q(i) + qDot(i) * dt;
            r(i+1) = r(i) + rDot(i) * dt;
            z(i+1) = z(i) + zdot(i) * dt;
            pitch(i+1) = pitch(i) + q(i) * dt;
            roll(i+1) = roll(i) + p(i) * dt;
            bank(i+1) = bank(i) + bdot(i) * dt;
            aoa(i+1) = aoa(i) + adot(i) * dt;
            ss(i+1) = ss(i) + ssdot(i) * dt;
            [lat(i+1), lon(i+i)] = propagateLatLon(lat(i), lon(i), V(i+1), heading(i+1), dt, constants);
            rotNedFromBody(:,:,i+1)  = angle2dcm(heading(i+1),pitch(i+1),roll(i+1),'ZYX');
            rotNedFromECEF(:,:,i+1)  = calcRotNedFromEcef(lat(i+1),lon(i+1));
            rotBodyFromECEF(:,:,i+1) = rotNedFromBody(:,:,i+1)' * rotNedFromECEF(:,:,i+1);
            chi(i+1) = calculateTrackAngle(V(i+1), lat(i+1), lon(i+1), heading(i+1), ss(i+1), constants, dt);

            % Simlog
            if i == (npoints-1)/4 || i == (npoints-1)/2 || i == 3 * (npoints-1)/4
                perc = 100 * i/(npoints-1);
                string1 = ['Simulation is ', round(num2str(perc)), '% complete'];
                disp(string1)
                alt
            end

            % End simulation when altitude reaches 0
            if alt <= 0
                break
            end
        end % end time loop
        dataSave.p(:, j, method + 1) = rad2deg(p);
        dataSave.q(:, j, method + 1) = rad2deg(q);
        dataSave.r(:, j, method + 1) = rad2deg(r);
        dataSave.bank(:, j, method + 1) = rad2deg(bank);
        dataSave.aoa(:, j, method + 1) = rad2deg(aoa);
        dataSave.heading(:, j, method + 1) = rad2deg(heading);
        dataSave.flightPathAngle(:, j, method + 1) = rad2deg(flightPathAngle);
        dataSave.gammaDot(:, j, method + 1) = rad2deg(gammaDot);
        dataSave.flightPathCmd(:, j, method + 1) = rad2deg(cmd_g(1,:));
        dataSave.headingCmd(:, j, method + 1) = rad2deg(cmd_g(2,:));
        dataSave.bankCmd(:, j, method + 1) = rad2deg(cmd_a(1,1:length(aoa)));
        dataSave.aoaCmd(:, j, method + 1) = rad2deg(cmd_a(2,1:length(aoa)));
        dataSave.ssCmd(:, j, method + 1) = rad2deg(cmd_a(3,1:length(aoa)));
        dataSave.pCmd(:, j, method + 1) = rad2deg(cmd_r(1,:));
        dataSave.qCmd(:, j, method + 1) = rad2deg(cmd_r(2,:));
        dataSave.rCmd(:, j, method + 1) = rad2deg(cmd_r(3,:));
        dataSave.pDotCmd(:, j, method + 1) = rad2deg(flin(1,:));
        dataSave.qDotCmd(:, j, method + 1) = rad2deg(flin(2,:));
        dataSave.rDotCmd(:, j, method + 1) = rad2deg(flin(3,:));
        dataSave.sideslip(:, method + 1) = rad2deg(ss);
        dataSave.time(:,j, method + 1) = tvec;
        fdplot = [flapDeflection(2,:);flapDeflection(3,:);flapDeflection(4,:);flapDeflection(5,:);flapDeflection(6,:);flapDeflection(7,:);flapDeflection(8,:);flapDeflection(1,:)];
        dataSave.fin1(:, j, method + 1) = fdplot(1,:);
        dataSave.fin2(:, j, method + 1) = fdplot(2,:);
        dataSave.fin3(:, j, method + 1) = fdplot(3,:);
        dataSave.fin4(:, j, method + 1) = fdplot(4,:);
        dataSave.fin5(:, j, method + 1) = fdplot(5,:);
        dataSave.fin6(:, j, method + 1) = fdplot(6,:);
        dataSave.fin7(:, j, method + 1) = fdplot(7,:);
        dataSave.fin8(:, j, method + 1) = fdplot(8,:);
        dataSave.fins(:,:,j, method + 1) = fdplot;

        disp('Simulation Complete')
        toc

        savestring = ['loopsel_' num2str(ctrl.loopselect) '_method' num2str(method)...
            'M' num2str(M) '_alt' num2str(alt_init/1000)];
        save(savestring)

        %% Plots
        set(0, 'DefaultLineLineWidth', 2);
        set(groot,'defaultAxesXGrid','on')
        set(groot,'defaultAxesYGrid','on')
    end % end method loop

%         if ctrl.loopselect == 0
%             figure('Name', ['Flight Path Angle and Heading' num2str(gyroNoiseDensity(j))]);
%             subplot(2,1,1); plot(tvec(2:end-1),180/pi * cmd_g(1,2:end),tvec,flightPathAngle * 180/pi);
%             ylabel('Flight Path Angle (deg)')
%             legend('Guidance Command','Response')
%             subplot(2,1,2); plot(tvec(2:end-1),180/pi * cmd_g(2,2:end),tvec,heading * 180/pi)
%             ylabel('headinging Angle (deg)'); xlabel('Time(sec)'); grid on
%
%             figure('Name', ['Altitude and Velocity' num2str(gyroNoiseDensity(j))])
%             subplot(2,1,1); plot(tvec,V)
%             title(['Flight Conditions with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('Velocity (m/s)')
%             subplot(2,1,2); plot(tvec,z-constants.RADIUS_EQ)
%             ylabel('Altitude (m)');xlabel('Time(sec)'); grid on
%         end
%
%         if ctrl.loopselect <= 1
%             figure('Name', ['Angle Outputs' num2str(gyroNoiseDensity(j))]);
%             subplot(3,1,1); plot(tvec,180/pi * cmd_a(1,1:length(aoa)),tvec,bank * 180/pi)
%             title(['Angle Output with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('bank (deg)')
%             legend('Angle Command','Response')
%             subplot(3,1,2); plot(tvec,180/pi * cmd_a(2,1:length(aoa)),tvec,aoa * 180/pi);
%             ylabel('aoa (deg)')
%             subplot(3,1,3); plot(tvec,180/pi * cmd_a(3,1:length(aoa)),tvec,ss * 180/pi)
%             ylabel('ss (deg)'); xlabel('Time(sec)'); grid on
%         end
%
%         figure('Name', ['Angular Rate Outputs' num2str(gyroNoiseDensity(j))]);
%         subplot(3,1,1); plot(tvec(1:end-1),cmd_r(1,:) * 180/pi,tvec,p * 180/pi)
%         title(['Angular Rate Output with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('p (deg/sec)')
%         legend('Command','Response')
%         subplot(3,1,2); plot(tvec(1:end-1),cmd_r(2,:) * 180/pi,tvec,q * 180/pi);
%         ylabel('q (deg/sec)')
%         subplot(3,1,3); plot(tvec(1:end-1),cmd_r(3,:) * 180/pi,tvec,r * 180/pi)
%         ylabel('r (deg/sec)'); xlabel('Time(sec)'); grid on
%
%         % figure;
%         % subplot(3,1,1); plot(tvec,p * 180/pi)
%         % title('Angular Rate Output'); ylabel('p (deg/sec)')
%         % subplot(3,1,2); plot(tvec,q * 180/pi);
%         % ylabel('q (deg/sec)')
%         % subplot(3,1,3); plot(tvec,r * 180/pi)
%         % ylabel('r (deg/sec)'); xlabel('Time(sec)'); grid on
%
%         % figure;
%         % subplot(3,1,1); plot(tvec,bank * 180/pi)
%         % title('Angle Displacements'); ylabel('Bank Angle (deg)')
%         % subplot(3,1,2); plot(tvec,aoa * 180/pi);
%         % ylabel('Angle of Attack (deg)')
%         % subplot(3,1,3); plot(tvec,ss * 180/pi)
%         % ylabel('Sideslip Angle (deg)'); xlabel('Time(sec)'); grid on
%
%         figure('Name', ['Angular Accel Output' num2str(gyroNoiseDensity(j))]);
%         subplot(3,1,1); plot(tvec(1:end-1),flin(1,:) * 180/pi,tvec(1:end-1),pDot(1:end-1) * 180/pi)
%         title(['Angular Accel Output with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('pDot (deg/s^2)')
%         % ylim([-2 2])
%         legend('Accel Command','Response')
%         subplot(3,1,2); plot(tvec(1:end-1),flin(2,:) * 180/pi,tvec(1:end-1),qDot(1:end-1) * 180/pi);
%         ylabel('qDot (deg/s^2)')
%         % ylim([-2 2])
%         subplot(3,1,3); plot(tvec(1:end-1),flin(3,:) * 180/pi, tvec(1:end-1),rDot(1:end-1) * 180/pi)
%         ylabel('rDot (deg/s^2)'); xlabel('Time(sec)'); grid on
%         % ylim([-2 2])
%
%         % figure;
%         % subplot(3,1,1); plot(tvec(1:end-1),pDot(1:end-1) * 180/pi)
%         % title('Angular Accel Output'); ylabel('pDot (deg/s^2)')
%         % subplot(3,1,2); plot(tvec(1:end-1),qDot(1:end-1) * 180/pi);
%         % ylabel('qDot (deg/s^2)')
%         % subplot(3,1,3); plot(tvec(1:end-1),rDot(1:end-1) * 180/pi)
%         % ylabel('rDot (deg/s^2)'); xlabel('Time(sec)'); grid on
%
%         % reorder flaps to match report numbering
%         fdplot = [flapDeflection(2,:);flapDeflection(3,:);flapDeflection(4,:);flapDeflection(5,:);flapDeflection(6,:);flapDeflection(7,:);flapDeflection(8,:);flapDeflection(1,:)];
%         dataSave.finDeflectSave(:,:,j) = fdplot;
%
%         figure('Name', ['Flap Deflections ' num2str(gyroNoiseDensity(j))])
%         plot(tvec(1:end-1),fdplot(:,1:end-1))
%         title(['Flap Deflect Angles with Gyro noise of ' num2str(gyroNoiseDensity(j))]); ylabel('deg')
%         legend('F1','F2','F3','F4','F5','F6','F7','F8')
%         xlabel('Time (Sec)')
%         grid on

end % end of gyro noise loop


% plotGyroNoise(dataSave)

% figure
% plot(tvec,gammaDot,tvec,gammaDotFilt)
% figure
% plot(tvec,headingDot,tvec,headingDotFilt)
%
% figure
% plot(tvec,pDot,tvec,pDotFilt)
% figure
% plot(tvec,qDot,tvec,qDotFilt)
% figure
% plot(tvec,rDot,tvec,rDotFilt)

% figure
% subplot(3,1,1)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))
% legend('raw','filtered')
% subplot(3,1,2)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))
% subplot(3,1,3)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))

%% Plots
controllers = {'NDI', 'INDI'};
colors = lines(length(controllers) + 2 * length(gyroNoiseDensity) - 1);
linS = {'-', '--'};
colors(end,:) = 255 / 255 * ones(1,3);
mode = 'Rate';
%% Angle Plots
figure('Name', "Angles")
tiledlayout(3,1)
nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(:,j,i), dataSave.bank(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Bank Angle (deg)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(:,j,i),dataSave.bankCmd(:,1,1), 'color', colors(count,:), 'linew', 2)
Legend{count} = 'Cmd';
legend(Legend, 'location', 'northoutside')

nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(:,j,i), dataSave.aoa(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Angle of Attack (deg)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(:,j,i), dataSave.aoaCmd(:,1,1), 'color', colors(count,:), 'linew', 2)
Legend{count} = 'AoA Cmd';
% legend(Legend, 'location', 'eastoutside')

nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(:,j,i), dataSave.heading(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Sideslip Angle (deg)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(:,1,1), dataSave.ssCmd(:,1,1), 'color', colors(count,:), 'linew', 2)
Legend{count} = 'SS Cmd';
% legend(Legend, 'location', 'eastoutside')w
% set(gcf, 'Position', get(0,'Screensize'));
if ~exist('/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/', 'dir')
    mkdir('/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/')
end
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/Angles' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/Angles' mode '.png'])

%% Flight Path Angle and Heading Plots
figure('Name', "Gamma and Heading")
tiledlayout(2,1)
nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(:,j,i), dataSave.flightPathAngle(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Flight Path Angle (deg)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(1:end-1,j,i),dataSave.flightPathCmd(:,1,1), 'color', colors(count,:), 'linew', 2)
Legend{count} = 'Gamma Cmd';
legend(Legend, 'location', 'northoutside')

nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(:,j,i), dataSave.heading(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Heading (deg)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(1:end-1,j,i), dataSave.headingCmd(:,1,1), 'color', colors(count,:), 'linew', 2)
Legend{count} = 'Heading Cmd';
% legend(Legend, 'location', 'eastoutside')

%% Angular Rate Plots
figure('Name', "Angle Rates")
tiledlayout(3,1)
nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(:,j,i), dataSave.p(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Roll Rate (deg/s)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(1:end-1,j,i), dataSave.pCmd(:,j,i), 'linew', 2, 'color', colors(count,:))
Legend{count} = 'Cmd';
legend(Legend, 'location', 'northoutside')

nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(:,j,i), dataSave.q(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Pitch Rate (deg/s)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(1:end-1,j,i), dataSave.qCmd(:,j,i), 'linew', 2, 'color', colors(count,:))
Legend{count} = 'Pitch Rate Cmd';
% legend(Legend, 'location', 'eastoutside')

nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(:,j,i), dataSave.r(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Yaw Rate (deg/s)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(1:end-1,j,i), dataSave.rCmd(:,j,i), 'linew', 2, 'color', colors(count,:))
Legend{count} = 'Yaw Rate Cmd';
% legend(Legend, 'location', 'eastoutside')
% set(gcf, 'Position', get(0,'Screensize'));
if ~exist('/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/')
    mkdir('/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/')
end
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/Rates' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/Rates' mode '.png'])

%% Angular Accel
figure('Name', "Angle Accel")
tiledlayout(3,1)
nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(1:end-1,j,i), dataSave.pDotSave(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Roll Accel (deg/s)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(1:end-1,j,i), dataSave.pDotCmd(:,j,i), 'linew', 2, 'color', colors(count,:))
Legend{count} = 'Roll Accel Cmd';
legend(Legend, 'location', 'eastoutside')

nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(1:end-1,j,i), dataSave.qDotSave(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Pitch Accel (deg/s)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(1:end-1,j,i), dataSave.qDotCmd(:,j,i), 'linew', 2, 'color', colors(count,:))
Legend{count} = 'Pitch Accel Cmd';
legend(Legend, 'location', 'eastoutside')

nexttile
hold on
count = 1;
for i = 1:length(methods)
    for j = 1:length(gyroNoiseDensity)
        plot(dataSave.time(1:end-1,j,i), dataSave.rDotSave(:,j,i), 'linew', 2, 'color', colors(count,:), 'linestyle', linS{i})
        Legend{count} = ['Noise density ' num2str(gyroNoiseDensity(j)) ' ' controllers{i}];
        ylabel('Yaw Accel (deg/s)')
        xlabel('Time (s)')
        count = count + 1;
    end
end
plot(dataSave.time(1:end-1,j,i), dataSave.rDotCmd(:,j,i), 'linew', 2, 'color', colors(count,:))
Legend{count} = 'Yaw Accel Cmd';
legend(Legend, 'location', 'eastoutside')
% set(gcf, 'Position', get(0,'Screensize'));
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/AngularAccel' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/AngularAccel' mode '.png'])

%% Fin Cmds
figure
plot(dataSave.time(:,1,1) .*  ones(size(dataSave.fins(:,:,1,1)')), dataSave.fins(:,:,1,1)', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Fin Cmds (deg)')
for i = 1:8
    Legend{i} = ['Fin' num2str(i)];
end
legend(Legend, 'location', 'eastoutside')
% set(gcf, 'Position', get(0,'Screensize'));
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/FinsNoNoiseNDI' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/FinsNoNoiseNDI' mode '.png'])

figure
plot(dataSave.time(:,1,1) .*  ones(size(dataSave.fins(:,:,1,1)')),dataSave.fins(:,:,1,2)', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Fin Cmds (deg)')
for i = 1:8
    Legend{i} = ['Fin' num2str(i)];
end
legend(Legend, 'location', 'eastoutside')
% set(gcf, 'Position', get(0,'Screensize'));
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/FinsNoNoiseINDI' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/FinsNoNoiseINDI' mode '.png'])

figure
plot(dataSave.time(:,1,1) .*  ones(size(dataSave.fins(:,:,1,1)')),dataSave.fins(:,:,end,1)', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Fin Cmds (deg)')
for i = 1:8
    Legend{i} = ['Fin' num2str(i)];
end
legend(Legend, 'location', 'eastoutside')
% set(gcf, 'Position', get(0,'Screensize'));
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/FinsMaxNoiseNDI' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/FinsMaxNoiseNDI' mode '.png'])

figure
plot(dataSave.time(:,1,1) .*  ones(size(dataSave.fins(:,:,1,1)')), dataSave.fins(:,:,end,2)', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Fin Cmds (deg)')
for i = 1:8
    Legend{i} = ['Fin' num2str(i)];
end
legend(Legend, 'location', 'eastoutside')
% set(gcf, 'Position', get(0,'Screensize'));
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/FinsMaxNoiseINDI' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/FinsMaxNoiseINDI' mode '.png'])

figure
plot(dataSave.time(:,1,1) .*  ones(size(dataSave.fins(:,:,1,1)')),dataSave.fins(:,:,end-1,1)', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Fin Cmds (deg)')
for i = 1:8
    Legend{i} = ['Fin' num2str(i)];
end
legend(Legend, 'location', 'eastoutside')
% set(gcf, 'Position', get(0,'Screensize'));
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/Fins1e4NoiseNDI' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/Fins1e4NoiseNDI' mode '.png'])

figure
plot(dataSave.time(:,1,1) .*  ones(size(dataSave.fins(:,:,1,1)')), dataSave.fins(:,:,end-1,2)', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('Fin Cmds (deg)')
for i = 1:8
    Legend{i} = ['Fin' num2str(i)];
end
legend(Legend, 'location', 'eastoutside')
% set(gcf, 'Position', get(0,'Screensize'));
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/Fins1e4NoiseINDI' mode '.fig'])
export_fig(['/home/amehdipour/repos/MastersProject/Sim/Plots/GammaSweep/Fins1e4NoiseINDI' mode '.png'])
