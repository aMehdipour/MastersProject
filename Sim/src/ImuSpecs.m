%==========================================================================
% This code takes specs from various navigation grade IMUs used in
% spacecraft to be used in a noise sensitivity analysis in the ADEPT
% simulation
%==========================================================================

%==========================================================================
% Honeywell HG4930:
%   Gyro bias stability: ? 0.002°/hr
%   Gyro noise density: ? 0.0006°/?hr or 3.053e-9 rad/s/?Hz
%   Accelerometer bias stability: ? 2 µg
%   Accelerometer noise density: ? 70 µg/?Hz
%==========================================================================

imu.HG4930.gyroBiasStab = rad2deg(0.002) / 3600;
imu.HG4930.gyroNoiseDensity = 3.053e-9;

%==========================================================================
% Northrop Grumman LN-200:
%   Gyro bias stability: ? 0.003°/hr
%   Gyro noise density: ? 0.003°/?hr
%   Accelerometer bias stability: ? 2 mg
%   Accelerometer noise density: ? 200 µg/?Hz
%==========================================================================

imu.LN200.gyroBiasStab = rad2deg(0.003)

%==========================================================================
% KVH Industries DSP-1760:
%   Gyro bias stability: ? 0.05°/hr
%   Gyro noise density: ? 0.004°/?hr
%   Accelerometer bias stability: ? 10 µg
%   Accelerometer noise density: ? 180 µg/?Hz
%==========================================================================

