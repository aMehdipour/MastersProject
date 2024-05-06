function [crossrange] = calculateCrossrange(state, constants)
    [~, terminalAzimuth] = distance('gc', rad2deg(state.latitude), rad2deg(state.longitude), constants.TARGET_LAT, constants.TARGET_LON);
    % lat = rad2deg(state.latitude);
    % lon = rad2deg(state.longitude);
    % X = cosd(constants.TARGET_LAT) * sind(constants.TARGET_LON - lon);
    % Y = cosd(lon) * sind(constants.TARGET_LON) - sind(lon) * cosd(constants.TARGET_LON) * cosd(constants.TARGET_LON - lon);
    % terminalAzimuth = atan2d(X,Y);

    terminalAzimuth = deg2rad(terminalAzimuth);
    crossrange = asin(sin(state.rangeToGo) * sin(state.heading - terminalAzimuth));
