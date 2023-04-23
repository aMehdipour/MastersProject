function [L,D,S,Lm,Mm,Nm] = ComputeAero(V,rho,aoa,ss,bank,fd,CFD,constants)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the aerodynamic loading on the entry vehicle based
% on the provided flight conditions and CFD data.
%
% Inputs
% V       vehicle velocity (m/s)
% rho     atmospheric density (kg/m^3)
% aoa     angle of attack (deg)
% ss      sideslip angle (deg)
% bank    bank angle (deg)
% fd      flap deflection angles arranged in a 1x8 vector (deg)
% CFD     CFD data structure
%
% Outputs
% L,D,S     Lift, drag and side forces (N)
% Lm,Mm,Nm  Roll, pitch and yaw moments (Nm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

err = CFD.err;

%% Compute vehicle forces and moments
ss   = ss * 180/pi;
aoa  = aoa * 180/pi;
avec = CFD.aoa(:,1); % vector of angle of attack datapoints
svec = CFD.ss(1,:);  % vector of sideslip angle datapoints
F_scale_body = (1/2)*rho * constants.FRONTAL_AREA * V^2; % scaling factor for forces (N)
M_scale_body = (1/2)*rho * constants.FRONTAL_AREA * constants.REFERENCE_LENGTH * V^2; % scaling factor for moments (Nm)

L_A      = F_scale_body * interp2(avec, svec, CFD.A.Cl,ss,aoa) * err; % query points ordering is flipped
D_A      = F_scale_body * interp2(avec, svec, CFD.A.Cd,ss,aoa) * err; % FIXME: Originally was not multiplied by err
S_A      = F_scale_body * interp2(avec, svec, CFD.A.Cs,ss,aoa) * err;
Lm_A_pre = M_scale_body * interp2(avec, svec, CFD.A.CmL,ss,aoa) * err; % aero moments computed in aerodynamic frame
Mm_A_pre = M_scale_body * interp2(avec, svec, CFD.A.CmM,ss,aoa) * err;
Nm_A_pre = M_scale_body * interp2(avec, svec, CFD.A.CmN,ss,aoa) * err;

% Convert moments to body frame through the bank angle
rotBodyFromBank = [1 0 0; 0 cos(bank) -sin(bank); 0 sin(bank) cos(bank)];
tmp  = rotBodyFromBank * [Lm_A_pre;Mm_A_pre;Nm_A_pre];
Lm_A = tmp(1);
Mm_A = tmp(2);
Nm_A = tmp(3);

%% Compute flap forces and moments
F_scale_flap = (1/2) * rho * constants.FLAP_AREA * V^2;
M_scale_flap = (1/2) * rho * constants.FLAP_AREA * constants.FLAP_MOMENT_ARM*V^2;

% initialize loop variables as zeros
L_F = zeros(1,8); D_F=L_F; S_F=L_F; Lm_F=L_F; Mm_F=L_F; Nm_F=L_F;

for i = 1:8

    if abs(fd(i)) > 20
        fd(i) = 20*sign(fd(i));
    end

    % set up grid points for data entries
    if i == 1
    X = CFD.ss3d; Y = CFD.aoa3d; Z = CFD.deflect3d;
    end

    ind = ['CFD.F' num2str(i) '.']; % flap indicator
    L_F(i)      = F_scale_flap * interp3(X, Y, Z, eval([ind 'Clline']),ss,aoa,fd(i)) * err;
    D_F(i)      = F_scale_flap * interp3(X, Y, Z, eval([ind 'Cdline']),ss,aoa,fd(i)) * err;
    S_F(i)      = F_scale_flap * interp3(X, Y, Z, eval([ind 'Csline']),ss,aoa,fd(i)) * err;
    Lm_F_pre(i) = M_scale_flap * interp3(X, Y, Z, eval([ind 'CmLline']),ss,aoa,fd(i)) * err;
    Mm_F_pre(i) = M_scale_flap * interp3(X, Y, Z, eval([ind 'CmMline']),ss,aoa,fd(i)) * err;
    Nm_F_pre(i) = M_scale_flap * interp3(X, Y, Z, eval([ind 'CmNline']),ss,aoa,fd(i)) * err;

    % Convert moments to body frame
    tmp = rotBodyFromBank * [Lm_F_pre(i);Mm_F_pre(i); Nm_F_pre(i)];
    Lm_F(i) = tmp(1);
    Mm_F(i) = tmp(2);
    Nm_F(i) = tmp(3);
end

%% Summation of Forces and Moments
L  = L_A + sum(L_F); % lift force (N)
D  = D_A + sum(D_F); % drag force (N)
S  = S_A + sum(S_F); % side force (N)
Lm = Lm_A + sum(Lm_F); % roll moment (Nm)
Mm = Mm_A + sum(Mm_F); % pitch moment (Nm)
Nm = Nm_A + sum(Nm_F); % yaw moment (Nm)
