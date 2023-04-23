function [v_aero,B,bias,f_aero,Ffit] = ControllerAero(V,rho,aoa,ss,bank,fd,CFD,constants)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the aerodynamic inputs to the control system based
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
% v_aero    Aerodynamic moments acting solely on the vehicle
% f_aero    Aerodynamic forces acting solely on the control surfaces
% B         Control effectiveness matrix, a 5x8 matrix containing the
%           slopes of the moment data with respect to flap deflection angles
% bias      The flap aerodynamic data is approximated using a linear
%           regression model. The slope of the linear regression is used to
%           construct the B matrix, and the bias is subtracted from
%           the final moment command to correct for the offset.
% Ffit      Linear fit (bias and slope) of vehicle forces with respect to
%            changes in angle of attack and sideslip angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute vehicle moments
ss = ss * 180/pi; aoa = aoa * 180/pi; % convert angles to degrees
aoaVec = CFD.aoa(:,1); % vector of angle of attack datapoints
ssVec = CFD.ss(1,:); % vector of sideslip angle datapoints
F_const = (1/2) * rho * constants.FRONTAL_AREA * V^2; % scaling factor for forces (N)
M_const = (1/2) * rho * constants.FRONTAL_AREA * constants.REFERENCE_LENGTH * V^2; % scaling factor for moments (Nm)

Lm_A_pre = M_const * interp2(aoaVec, ssVec, CFD.A.CmL, ss, aoa); %aero moments computed in aerodynamic frame
Mm_A_pre = M_const * interp2(aoaVec, ssVec, CFD.A.CmM, ss, aoa);
Nm_A_pre = M_const * interp2(aoaVec, ssVec, CFD.A.CmN, ss, aoa);

% Convert moments to body frame through the bank angle
rotBodyFromBank = [1 0 0; 0 cos(bank) -sin(bank); 0 sin(bank) cos(bank)];
tmp  = rotBodyFromBank * [Lm_A_pre;Mm_A_pre;Nm_A_pre];
Lm_A = tmp(1);
Mm_A = tmp(2);
Nm_A = tmp(3);

% Store vehicle moments for NDI
v_aero = struct('Lmv',Lm_A,'Mmv',Mm_A,'Nmv',Nm_A);

% Fit lift and side force data to range of angle deltas
angvec = [-20 -10 0 10 20];
for i = 1:5
    CL_A(i) = interp2(aoaVec, ssVec, CFD.A.Cl,ss,angvec(i)); %query points ordering is flipped
    CS_A(i) = interp2(aoaVec, ssVec, CFD.A.Cs,angvec(i),aoa);
end

Clfit = polyfit(-20:10:20, CL_A,1);
Csfit = polyfit(-20:10:20, CS_A,1);
Ffit  = F_const * [Clfit(2) Csfit(2) Clfit(1) Csfit(1)]; % [bias slope]

%% Compute flap lift and drag forces
F_const = (1/2) * rho * constants.FLAP_AREA * V^2;
M_const = (1/2) * rho * constants.FLAP_AREA * constants.FLAP_MOMENT_ARM * V^2; % scaling factor for moments (Nm)

% initialize loop variables as zeros
L_F = zeros(1,8); S_F=L_F;

for i = 1:8

    if abs(fd(i)) > 20
        fd(i) = 20 * sign(fd(i));
    end

    % set up grid points for data entries
    if i == 1
    X = CFD.ss3d; Y = CFD.aoa3d; Z = CFD.deflect3d;
    end

    ind = ['CFD.F' num2str(i) '.']; % flap indicator
    L_F(i) = F_const * interp3(X, Y, Z, eval([ind 'Clline']),ss,aoa,fd(i));
    S_F(i) = F_const * interp3(X, Y, Z, eval([ind 'Csline']),ss,aoa,fd(i));
end

L_f = sum(L_F); % lift force (N)
S_f = sum(S_F); % side force (N)

% Store flap forces for NDI
f_aero = struct('Lf',L_f,'Sf',S_f);

%% Interpolate slope and offset data for B matrix and bias computation
bias = zeros(3,1);
B = zeros(3,8);
for i = 1:8
    ind = ['CFD.F', num2str(i), '.'];
    CmLslope = M_const * interp2(aoaVec, ssVec,  eval([ind,'CmLslope']),ss,aoa);
    CmMslope = M_const * interp2(aoaVec, ssVec,  eval([ind,'CmMslope']),ss,aoa);
    CmNslope = M_const * interp2(aoaVec, ssVec,  eval([ind,'CmNslope']),ss,aoa);

    B(:,i) = rotBodyFromBank * [CmLslope; CmMslope; CmNslope];

    CmLbias = M_const * interp2(aoaVec, ssVec, eval([ind,'CmLbias']),ss,aoa);
    CmMbias = M_const * interp2(aoaVec, ssVec, eval([ind,'CmMbias']),ss,aoa);
    CmNbias = M_const * interp2(aoaVec, ssVec, eval([ind,'CmNbias']),ss,aoa);

    bias = bias + rotBodyFromBank * [CmLbias; CmMbias; CmNbias];
end
