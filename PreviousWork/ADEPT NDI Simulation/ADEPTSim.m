%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADEPTSim
% This script simulates the atmospheric entry of an ADEPT-class entry
% vehicle utilizing an aerodynamic flap control system.
%
% Created By: Joshua Stokes
% Last Modified: 8/14/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all
addpath Data; addpath('QCAT\qcat'); addpath src; addpath Outputs
load CFD; load std_atmos; load att_profile

%% Constants
c = 271.42; % speed of sound at 100km altitude (m/s)
Rearth = 6378000; % radius of Earth (m)
m = 8.490; % vehicle mass (kg)
I = [0.25 0 0; 0 0.1722 0; 0 0 0.1719]; % vehicle moment of inertia tensor (kg-m^2)

%% Sim Settings
dt = 1/500;
tend = 5;
tvec = 0:dt:tend;
npoints = length(tvec);
CFD.err = 1;

%% Controller
ctrl.ts = 1/500; % sample time (s)
ts_g = 1/500;
ctrl.tsratio = ts_g/ctrl.ts;
ctrl.cnt = 1;
method = 0; % 0 to use NDI,1 to use INDI
ctrl.loopselect = 1; % 0 to use guidance, 1 to use angle, 2 to use rate
cs_rate = 200; % maximum deflection rate for control surfaces (deg/s)
ctrl.maxdeflect = ctrl.ts*cs_rate; %maximum deflection over one ctrl sample (deg)
ctrl.ratelim = 10*pi/180;
ctrl.anglim = 18*pi/180;

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
OL_rate = ctrl.kp_r*DI;
CL_rate = feedback(OL_rate,1); % closed rate loop

OL_angle = ctrl.kp_a*CL_rate;
CL_angle = feedback(OL_angle,1);

% pidTuner(CL_rate)
% pidTuner(DI*CL_angle)
CNT_guidance = tf([ctrl.kd_g ctrl.kp_g ctrl.ki_g],[0 1 0]);

% Active set inputs
ctrl.Wv = [1 0 0; 0 1 0; 0 0 1];%eye(3);
ctrl.Wu = eye(8);
ctrl.gam_wls = 1e6;

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
gamma = zeros(npoints,1); gammadot = gamma; % flight path angle (rad)
head = zeros(npoints,1); hdot = head; % heading angle (rad)
p = zeros(npoints,1); pdot = p; % angular velocity about x axis (rad/s)
q = zeros(npoints,1); qdot = q; % angular velocity about y axis (rad/s)
r = zeros(npoints,1); rdot = r; % angular velocity about z axis (rad/s)
pdotfilt = r; qdotfilt = r; rdotfilt = r; 
gdotfilt = r; hdotfilt = r;

% Kinematics
z = zeros(npoints,1); zdot = z; % orbital radius (m)
lat = zeros(npoints,1); latdot = lat; % latitude (rad)
lon = zeros(npoints,1); londot = lon; % longitude (rad)
bank = zeros(npoints,1); bdot = bank; % bank angle (rad)
aoa = zeros(npoints,1); adot = aoa; % angle of attack (rad)
ss = zeros(npoints,1); ssdot = ss; % sideslip angle (rad)

% Telemetry
vlin = zeros(3,npoints-1);
flin = zeros(3,npoints-1);
if ~exist('cmd_a'); cmd_a = zeros(3,npoints); end
cmd_r = zeros(3,npoints-1);
cmd_g = zeros(2,npoints-1);

%% Initial Conditions
M = 30; % Mach number (dless) and speed of sound (m/s)
V(1) = M*c; gamma(1) = -5.5*pi/180; head(1) = 0;
p(1) = 0; q(1) = 0; r(1) = 0;
alt = 100000; %118000; 
alt_init = alt; % save initial altitude for filename
z(1) = alt + Rearth;
lat(1) = 0; lon(1) = 0; 
bank(1) = 0; aoa(1) = 0; ss(1) = 0;

fd = zeros(8,npoints); % flap deflection angles (deg)
fc =zeros(8,1);

% controller states
ctrl.elast_g = [0; 0];
ctrl.eint_g = [0; 0];
ctrl.fclast = fc;

%% Case Specific Overrides
% This section overrides the vehicle initial conditions to result in trim
% conditions at the start of the simulation

if 0 %M30 A100 
    fc = [-0.7761   -0.9669   -1.4826    1.1879    1.5912    2.6582    0.3138   -0.0417]';
elseif 0 %M30 A100 FPA-5.5
    aoa(1) = 2.395*pi/180; ss(1) = -0.400*pi/180;
    fc = [1.6441    0.1347   -1.6763   -1.3564   -1.2157    0.5829    1.3549    1.8832]';
end
%% Execute Computations
tic
for i = 1:npoints-1
    
    % Interpolate density and gravitational force
    alt = z(i)-Rearth;
    rho = interp1(std.alt,std.rho,alt,'linear','extrap'); % atmospheric density (kg/m^3)
    gz = m*interp1(std.alt,std.g,alt,'linear','extrap'); % gravitational force (N)
    
    % Execute Controller 
    if mod(i*dt,ctrl.ts) == 0 && i >1
        
        % Compute aerodynamic force/moment inputs to controller
        [v_aero,B,bias,f_aero,Ffit] = ControllerAero(V(i),rho,aoa(i),ss(i),bank(i),fd(:,i-1),CFD);
        
        % Save state variables into input structure
        st = struct('g',gamma(i),'V',V(i),'head',head(i),'p',p(i),'q',q(i),'r',r(i),'bank',bank(i),'rho',rho,'fd',fd(:,i-1),...
            'gdot',gdotfilt(i-1),'z',z(i-1),'hdot',hdotfilt(i-1),'pdot',pdotfilt(i-1),'qdot',qdotfilt(i-1),'rdot',rdotfilt(i-1),...
            'aoa',aoa(i),'ssl',ss(i),'gz',gz);

        % Commands (Input selected based on ctrl.loopselect)        
        ctrl.cmd_g = [-5.5;0]*pi/180; % p, q, r
%         cmd_g(1,1) = -5.5*pi/180;
%         ctrl.cmd_g = [cmd_g(1,i-1);0] + [0.75*V(i)/z(i)*cos(gamma(i));0]*dt;
%         cmd_slope = 0.02*pi/180; %rad/sec
%         if i*dt < 10
%             ctrl.cmd_g = [gamma(1); head(1)] - cmd_slope*[i*dt; 0];
%         else
%             ctrl.cmd_g = [cmd_g(1,i-1); head(1)] + cmd_slope*[0; i*dt-10];
%         end
%         if i == 2
%             ctrl.cmd_a = [0; aoa(2); ss(2)];
%         end
           
        if ctrl.loopselect == 1 
            if i*dt > 1 && i*dt < 3
                ctrl.cmd_a = [0;1; 0]*pi/180;
            elseif i*dt > 1
                ctrl.cmd_a = [0; 1; 1]*pi/180;
            else
                ctrl.cmd_a = [0; 0; 0]*pi/180;
            end
%               ctrl.cmd_a = cmd_a(:,i);
        end
        
        if i*dt > 0.5 && i*dt < 1.5
            ctrl.cmd_r = [0; 1;0]*pi/180;
        elseif i*dt > 0.5
            ctrl.cmd_r = [0; 1;1]*pi/180;
        else
            ctrl.cmd_r = [0; 0; 0];
        end
        
        % Execute control law
        [fc,ctrl] = RunController(st,ctrl,v_aero,f_aero,Ffit,B,bias,method);

        %store telemetry
        vlin(:,i) = ctrl.vlin;
        flin(:,i) = ctrl.flin;
        cmd_a(:,i) = ctrl.cmd_a;
        cmd_r(:,i) = ctrl.cmd_r;
        cmd_g(:,i) = ctrl.cmd_g;
    end
    
    % Update flap Positions
    fd(:,i) = fc;

    % Compute aerodynamic forces and moments
    [L,D,S,Lm,Mm,Nm] = ComputeAero(V(i),rho,aoa(i),ss(i),bank(i),fd(:,i),CFD);

    % Compute State Derivatives
    Vdot(i) = -(gz/m)*sin(gamma(i)) - (1/m)*D;
    gammadot(i) = (-gz/(m*V(i)) + (V(i)/z(i)))*cos(gamma(i)) + 1/(m*V(i))*(L*cos(bank(i)) - S*sin(bank(i)));
    hdot(i) = (1/(m*V(i)*cos(gamma(i))))*(L*sin(bank(i)) + S*cos(bank(i)));
    
    gyro = cross([p(i); q(i); r(i)],I*[p(i); q(i); r(i)]); %gyroscopic torque (Nm)
    Omegadot = inv(I)*([Lm; Mm; Nm;] - gyro); % rotational dynamics derivatives
    pdot(i) = Omegadot(1);
    qdot(i) = Omegadot(2);
    rdot(i) = Omegadot(3);

    zdot(i) = V(i)*sin(gamma(i));    
    bdot(i) = p(i)*cos(aoa(i))*sec(ss(i)) + r(i)*sin(aoa(i))*sec(ss(i));
    adot(i) = -p(i)*cos(aoa(i))*tan(ss(i)) + q(i) - r(i)*sin(aoa(i))*tan(ss(i));
    ssdot(i) = p(i)*sin(aoa(i)) - r(i)*cos(aoa(i));
    
    %Filter derivatives for INDI processing
    % Store filter inputs [in(n) in(n-1) in(n-2)]
    p_in = [pdot(i) p_in(1:2)]; q_in = [qdot(i) q_in(1:2)]; r_in = [rdot(i) r_in(1:2)];
    g_in = [gammadot(i) g_in(1:2)]; h_in = [hdot(i) h_in(1:2)];
    % Execute filter function
    pdotfilt(i) = filter_discrete(p_in, p_mem,a,b);
    qdotfilt(i) = filter_discrete(q_in, q_mem,a,b);
    rdotfilt(i) = filter_discrete(r_in, r_mem,a,b);
    gdotfilt(i) = filter_discrete(g_in, g_mem,a_g,b_g);
    hdotfilt(i) = filter_discrete(h_in, h_mem,a_g,b_g);

    % Store filter outputs into "memory" [out(n) out(n-1)]
    p_mem = [pdotfilt(i) p_mem(1)]; q_mem = [qdotfilt(i) q_mem(1)]; r_mem = [rdotfilt(i) r_mem(1)];
    g_mem = [gdotfilt(i) g_mem(1)]; h_mem = [hdotfilt(i) h_mem(1)];

    % Update State Variables
    V(i+1) = V(i) + Vdot(i)*dt;
    gamma(i+1) = gamma(i) + gammadot(i)*dt;
    head(i+1) = head(i) + hdot(i)*dt;
    p(i+1) = p(i) + pdot(i)*dt;
    q(i+1) = q(i) + qdot(i)*dt;
    r(i+1) = r(i) + rdot(i)*dt;
    z(i+1) = z(i) + zdot(i)*dt;
    bank(i+1) = bank(i) + bdot(i)*dt;
    aoa(i+1) = aoa(i) + adot(i)*dt;
    ss(i+1) = ss(i) + ssdot(i)*dt;

    % Simlog
    if i == (npoints-1)/4 || i == (npoints-1)/2 || i == 3*(npoints-1)/4
        perc = 100*i/(npoints-1);
        string1 = ['Simulation is ', round(num2str(perc)), '% complete'];
        disp(string1)
        alt
    end
    
    % End simulation when altitude reaches 0
    if alt <= 0
        break
    end
end
disp('Simulation Complete')
toc

savestring = ['loopsel_' num2str(ctrl.loopselect) '_method' num2str(method)...
              'M' num2str(M) '_alt' num2str(alt_init/1000)];
save(savestring)

%% Plots
set(0, 'DefaultLineLineWidth', 1.3);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

if ctrl.loopselect == 0
    figure; 
    subplot(2,1,1); plot(tvec(2:end-1),180/pi*cmd_g(1,2:end),tvec,gamma*180/pi);
    ylabel('Flight Path Angle (deg)')
    legend('Guidance Command','Response')
    subplot(2,1,2); plot(tvec(2:end-1),180/pi*cmd_g(2,2:end),tvec,head*180/pi)
    ylabel('Heading Angle (deg)'); xlabel('Time(sec)'); grid on
    
    figure
    subplot(2,1,1); plot(tvec,V)
    title('Flight Conditions'); ylabel('Velocity (m/s)')
    subplot(2,1,2); plot(tvec,z-Rearth)
    ylabel('Altitude (m)');xlabel('Time(sec)'); grid on
end

if ctrl.loopselect <= 1
    figure;
    subplot(3,1,1); plot(tvec,180/pi*cmd_a(1,1:length(aoa)),tvec,bank*180/pi)
    title('Angle Output'); ylabel('bank (deg)')
    legend('Angle Command','Response')
    subplot(3,1,2); plot(tvec,180/pi*cmd_a(2,1:length(aoa)),tvec,aoa*180/pi);
    ylabel('aoa (deg)')
    subplot(3,1,3); plot(tvec,180/pi*cmd_a(3,1:length(aoa)),tvec,ss*180/pi)
    ylabel('ss (deg)'); xlabel('Time(sec)'); grid on
end

figure;
subplot(3,1,1); plot(tvec(1:end-1),cmd_r(1,:)*180/pi,tvec,p*180/pi)
title('Angular Rate Output'); ylabel('p (deg/sec)')
legend('Command','Response')
subplot(3,1,2); plot(tvec(1:end-1),cmd_r(2,:)*180/pi,tvec,q*180/pi);
ylabel('q (deg/sec)')
subplot(3,1,3); plot(tvec(1:end-1),cmd_r(3,:)*180/pi,tvec,r*180/pi)
ylabel('r (deg/sec)'); xlabel('Time(sec)'); grid on

% figure;
% subplot(3,1,1); plot(tvec,p*180/pi)
% title('Angular Rate Output'); ylabel('p (deg/sec)')
% subplot(3,1,2); plot(tvec,q*180/pi);
% ylabel('q (deg/sec)')
% subplot(3,1,3); plot(tvec,r*180/pi)
% ylabel('r (deg/sec)'); xlabel('Time(sec)'); grid on

% figure;
% subplot(3,1,1); plot(tvec,bank*180/pi)
% title('Angle Displacements'); ylabel('Bank Angle (deg)')
% subplot(3,1,2); plot(tvec,aoa*180/pi); 
% ylabel('Angle of Attack (deg)')
% subplot(3,1,3); plot(tvec,ss*180/pi)
% ylabel('Sideslip Angle (deg)'); xlabel('Time(sec)'); grid on

figure;
subplot(3,1,1); plot(tvec(1:end-1),flin(1,:)*180/pi,tvec(1:end-1),pdot(1:end-1)*180/pi)
title('Angular Accel Output'); ylabel('pdot (deg/s^2)')
% ylim([-2 2])
legend('Accel Command','Response')
subplot(3,1,2); plot(tvec(1:end-1),flin(2,:)*180/pi,tvec(1:end-1),qdot(1:end-1)*180/pi); 
ylabel('qdot (deg/s^2)')
% ylim([-2 2])
subplot(3,1,3); plot(tvec(1:end-1),flin(3,:)*180/pi, tvec(1:end-1),rdot(1:end-1)*180/pi)
ylabel('rdot (deg/s^2)'); xlabel('Time(sec)'); grid on
% ylim([-2 2])

% figure;
% subplot(3,1,1); plot(tvec(1:end-1),pdot(1:end-1)*180/pi)
% title('Angular Accel Output'); ylabel('pdot (deg/s^2)')
% subplot(3,1,2); plot(tvec(1:end-1),qdot(1:end-1)*180/pi); 
% ylabel('qdot (deg/s^2)')
% subplot(3,1,3); plot(tvec(1:end-1),rdot(1:end-1)*180/pi)
% ylabel('rdot (deg/s^2)'); xlabel('Time(sec)'); grid on

% reorder flaps to match report numbering
fdplot = [fd(2,:);fd(3,:);fd(4,:);fd(5,:);fd(6,:);fd(7,:);fd(8,:);fd(1,:)];

figure
plot(tvec(1:end-1),fdplot(:,1:end-1))
title('Flap Deflect Angles'); ylabel('deg')
legend('F1','F2','F3','F4','F5','F6','F7','F8')
xlabel('Time (Sec)')
grid on

% figure
% plot(tvec,gammadot,tvec,gdotfilt)
% figure
% plot(tvec,hdot,tvec,hdotfilt)
% 
% figure
% plot(tvec,pdot,tvec,pdotfilt)
% figure
% plot(tvec,qdot,tvec,qdotfilt)
% figure
% plot(tvec,rdot,tvec,rdotfilt)

% figure
% subplot(3,1,1)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))
% legend('raw','filtered')
% subplot(3,1,2)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))
% subplot(3,1,3)
% plot(tvec(2:end),vlin(1,:),tvec(2:end),flin(1,:))