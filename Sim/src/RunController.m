function [fc,ctrl] = RunController(st,ctrl,v_aero,f_aero,Ffit,B,bias,method,constants)

%% Constants
Ixx = constants.INERTIA(1,1);
Iyy = constants.INERTIA(2,2);
Izz = constants.INERTIA(3,3);

%% Unpack Structures
% Unpack structure inputs
fieldNames = fieldnames(st);
for i = 1:length(fieldNames)
    eval(strcat(fieldNames{i}, '= st.', fieldNames{i},';'))
end

fieldNames = fieldnames(v_aero);
for i = 1:length(fieldNames)
    eval(strcat(fieldNames{i}, '= v_aero.', fieldNames{i},';'))
end

fieldNames = fieldnames(ctrl);
for i = 1:length(fieldNames)
    eval(strcat(fieldNames{i}, '= ctrl.', fieldNames{i},';'))
end

Lf = f_aero.Lf; Sf = f_aero.Sf;

%% Guidance
% Compute error terms and output required derivatives
errorGammaHeading      = cmd_g - [flightPathAngle; heading];
errorDotGammaHeading   = (errorGammaHeading - errorLastGammaHeading) / ts;
errorIntegGammaHeading = errorIntegGammaHeading + errorGammaHeading * ts;

v_g = kd_g * errorDotGammaHeading + kp_g * errorGammaHeading + ki_g * errorIntegGammaHeading;

% Filter the control signal
controlSignalGammaFilt   = filter_discrete([v_g(1) gmem(1,:)], gfmem(1,:),a_g,b_g);
controlSignalHeadingFilt = filter_discrete([v_g(2) gmem(2,:)], gfmem(2,:),a_g,b_g);

gflin = [controlSignalGammaFilt(end);controlSignalHeadingFilt(end)];

gammaReq = gflin(1); headingReq = gflin(2);

fbias  = Ffit(1:2);
fslope = Ffit(3:4) * 180/pi; % convert from N/deg to N/rad

% Dynamic inversion
if method == 0 % NDI
    tmp(1,1)  = constants.MASS * V * gammaReq + (gz - (constants.MASS * V^2 / z)) * cos(flightPathAngle) - Lf * cos(bank) + Sf * sin(bank);
    tmp(2,1)  = constants.MASS * V * cos(flightPathAngle) * headingReq - Lf * sin(bank) - Sf * cos(bank);
    Brot      = [cos(bank) -sin(bank); sin(bank) cos(bank)];
    forcereqs = inv(Brot) * tmp; % [Lcs;Scs]

    cmd_a = (forcereqs - fbias') ./ fslope';
    cmd_a = [0;cmd_a];

elseif method == 1 % INDI
    tmp(1,1)   = constants.MASS * V * (gammaReq - gammaDot);
    tmp(2,1)   = constants.MASS * V * cos(flightPathAngle) * (headingReq - headingDot);
    Brot       = [cos(bank) -sin(bank); sin(bank) cos(bank)];
    dforcereqs = inv(Brot) * tmp; % [Lcs;Scs]

    dcmd_a = (dforcereqs) ./ fslope';
    cmd_a = [0; aoa; ssl] + [0; dcmd_a];
end

% Angle Limiter
angvec = [bank; aoa; ssl];
for i = 1:3
    if abs(cmd_a(i)) > anglim
        cmd_a(i) = anglim * sign(cmd_a(i));

        % Integrator anti-windup logic
        if i == 2 || i == 3 %implement on pitch and yaw axes
            if abs(angvec(i)) >= anglim - 0.1 * pi/180 % if
                errorIntegGammaHeading(i-1) = 0;
            end
        end

    end
end

%% Angle Controller
% Compute error terms and output angular rate commands
if loopselect == 1 % If we are in angle control mode
 cmd_a = ctrl.cmd_a;
end

err_a = cmd_a - [bank; aoa; ssl];
cmd_r = kp_a * err_a;
% cmd_r(3) = -cmd_r(3); %negative yaw rate corresponds to positive ss angle

cmd_r(1) = (cmd_r(1) - r * sin(aoa) * sec(ssl))/(cos(aoa) * sec(ssl));
cmd_r(2) = cmd_r(2) + p * cos(aoa) * tan(ssl) + r * sin(aoa) * tan(ssl);
cmd_r(3) = (-cmd_r(3) + p * sin(aoa))/cos(aoa);


% Rate Limiter
for i = 1:3
    if abs(cmd_r(i)) > ratelim
        cmd_r(i) = ratelim * sign(cmd_r(i));
    end
end

%% Rotational Controller
% Compute error terms and output angular acceleration commands
if loopselect == 2
 cmd_r = ctrl.cmd_r;
end

err_r = cmd_r - [p; q; r];
vlin = kp_r * err_r;

% Filter the control signal
vpfilt = filter_discrete([vlin(1) vmem(1,:)], fmem(1,:),a,b);
vqfilt = filter_discrete([vlin(2) vmem(2,:)], fmem(2,:),a,b);
vrfilt = filter_discrete([vlin(3) vmem(3,:)], fmem(3,:),a,b);

flin = [vpfilt(end);vqfilt(end);vrfilt(end)];
% flin = [0;1;1] * pi/180;

% Output filtered signals to dynamic inversion

% pReq = 0; qReq = 0; rReq =0;
% Dynamic Inversion
if method == 0 % NDI
    pReq = vlin(1); qReq = vlin(2); rReq = vlin(3);
    Lmcs = Ixx * pReq + Izz * r * q - Iyy * q * r - Lmv; % ignoring off-axis inertia (all zeros)
    Mmcs = Iyy * qReq + Ixx * p * r - Izz * p * r - Mmv;
    Nmcs = Izz * rReq - Ixx * p * q + Iyy * p * q - Nmv;
    v = [Lmcs; Mmcs; Nmcs] - bias;
elseif method == 1 % INDI
    pReq = flin(1); qReq = flin(2); rReq = flin(3);
    wReq = [pReq; qReq; rReq];
    wDot = [pDot; qDot; rDot];
    tmp2 = constants.INERTIA * (wReq - wDot);
    Lmcs = tmp2(1);
    Mmcs = tmp2(2);
    Nmcs = tmp2(3);
    v = [Lmcs; Mmcs; Nmcs];
end


% %% Control Allocation
% umin = zeros(8,1); dumin = zeros(8,1);
% umax = zeros(8,1); dumax = zeros(8,1);
% if method == 0
%     ud = zeros(8,1);

%     for i = 1:8
%         umin(i) = max([-constants.FLAP_LIMIT fclast(i) - maxdeflect]);
%         umax(i) = min([constants.FLAP_LIMIT fclast(i) + maxdeflect]);
% %         umin(i) = -20;
% %         umax(i) = 20;
%         if ud(i) < umin(i)
%             ud(i) = umin(i);
%         elseif ud(i) > umax(i)
%             ud(i) = umax(i);
%         end
%     end
%     % implement weighted least squares
%     [fc,W,iter] = wls_alloc(B,v,umin,umax,Wv,Wu,ud,gam_wls);

% elseif method == 1
%     dud = zeros(8,1) - flapDeflection;

%     for i = 1:8
%         dumin(i) = max([-constants.FLAP_LIMIT-flapDeflection(i) -maxdeflect]);
%         dumax(i) = min([constants.FLAP_LIMIT-flapDeflection(i) maxdeflect]);
%         if dud(i) < dumin(i)
%             dud(i) = dumin(i);
%         elseif dud(i) > dumax(i)
%             dud(i) = dumax(i);
%         end
%     end
%     % implement weighted least squares
%     [dfc,W,iter] = wls_alloc(B,v,dumin,dumax,Wv,Wu,dud);
%     fc = fclast + dfc;
% end

%% Control Allocation
umin = zeros(8,1); dumin = zeros(8,1);
umax = zeros(8,1); dumax = zeros(8,1);
if method == 0
    ud = zeros(8,1);

    for i = 1:8
        umin(i) = max([-constants.FLAP_LIMIT fclast(i) - maxdeflect]);
        umax(i) = min([constants.FLAP_LIMIT fclast(i) + maxdeflect]);
        if ud(i) < umin(i)
            ud(i) = umin(i);
        elseif ud(i) > umax(i)
            ud(i) = umax(i);
        end
    end
    % implement weighted least squares with PCH
    [fc,W,iter] = wls_alloc_pch(B,v,umin,umax,Wv,Wu,ud,gam_wls,(umin+umax)/2, zeros(length(umin),1), 100, 1e-6, 5);

elseif method == 1
    dud = zeros(8,1) - flapDeflection;

    for i = 1:8
        dumin(i) = max([-constants.FLAP_LIMIT-flapDeflection(i) -maxdeflect]);
        dumax(i) = min([constants.FLAP_LIMIT-flapDeflection(i) maxdeflect]);
        if dud(i) < dumin(i)
            dud(i) = dumin(i);
        elseif dud(i) > dumax(i)
            dud(i) = dumax(i);
        end
    end
    % implement weighted least squares with PCH
    [dfc,W,iter] = wls_alloc_pch(B,v,dumin,dumax,Wv,Wu,dud,gam_wls,(dumin+dumax)/2, zeros(length(dumin),1), 100, 1e-6, 5);
    fc = fclast + dfc;
end

%% Save Controller Data
ctrl.errorLastGammaHeading = errorGammaHeading;
ctrl.errorIntegGammaHeading = errorIntegGammaHeading;
ctrl.elast_a = err_a; % FIXME: Was err_r, was this correct? Changed to err_a
ctrl.elast_r = err_r;
ctrl.fclast = fc;
ctrl.vlin = vlin;
ctrl.flin = flin;
ctrl.vmem = [vlin ctrl.vmem(:,2)];
ctrl.fmem = [flin ctrl.fmem(:,1)];
ctrl.gmem = [v_g ctrl.gmem(:,2)];
ctrl.gfmem = [gflin ctrl.gfmem(:,1)];
ctrl.cmd_r = cmd_r;
ctrl.cmd_a = cmd_a;
ctrl.cnt = cnt;
