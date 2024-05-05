function headingError = calculateHeadingError(state, constants)
    % Calculate the azimuth between the current location and the target
    [~, azimuth] = distance('gc', state.latitude, state.longitude, constants.TARGET_LAT, constants.TARGET_LON);

    % Calculate the heading error
    headingError = state.heading - azimuth;
end%
