function [crossrange] = calculateCrossrange(state, constants)
    [~, terminalAzimuth] = distance('gc', state.latitude, state.longitude, constants.TARGET_LAT, constants.TARGET_LON);
    crossrange = asin(sin(state.rangeToGo * sin(state.heading - terminalAzimuth)));
