function [downrange, crossrange] = calculateDownrangeCrossrange(state, constants)
    % Constants
    deg2rad = pi / 180;

    % Convert target (trajectory system origin) latitude and longitude to ECEF coordinates
    tgtPosECEF = lla2ecef([constants.TARGET_LAT, constants.TARGET_LON, 0], 'WGS84');

    % Convert vehicle's current latitude, longitude, and altitude to ECEF coordinates
    vehiclePosECEF = lla2ecef([state.latitude, state.longitude, state.altitudeUnnormalized - constants.EARTH_RADIUS_EQ], 'WGS84');

    % Calculate ECEF position vector from target to vehicle
    posVehicleRelTgtECEF = vehiclePosECEF - tgtPosECEF;

    % Convert angles from degrees to radians
    targetLatRad = constants.TARGET_LAT * deg2rad;
    targetLonRad = constants.TARGET_LON * deg2rad;
    [~, azimuthRad] = distance("gc",rad2deg(state.latitude), rad2deg(state.longitude), constants.TARGET_LAT, constants.TARGET_LON); % Set the azimuth of the trajectory coordinate system, to be defined as needed
    azimuthRad = wrapTo2Pi(azimuthRad * deg2rad);

    % Step 2: Rotate by -targetLonRad about the W-axis to align the U-axis toward the polar axis
    R1 = [-cos(targetLonRad), -sin(targetLonRad), 0; -sin(targetLonRad), cos(targetLonRad), 0; 0, 0, 1];
    intermediate = R1 * posVehicleRelTgtECEF';

    % Step 3: Rotate by -targetLatRad about the U-axis to align the n-axis with the North-East (NE) plane
    R2 = [sin(targetLatRad), 0, cos(targetLatRad); -cos(targetLatRad), 0, sin(targetLatRad); 0, 1, 0];
    intermediate = R2 * intermediate;

    % Step 4: Rotate by -azimuthRad about the n-axis to align with the downrange-vertical-crossrange coordinate system
    R3 = [cos(azimuthRad), 0, sin(azimuthRad); 0, 1, 0; -sin(azimuthRad), 0, cos(azimuthRad)];
    final = R3 * intermediate;

    % Extract downrange and crossrange positions
    downrange = final(1);
    crossrange = final(3);
end
