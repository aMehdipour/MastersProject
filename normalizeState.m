function stateOut = normalizeState(state, normalizationType, flagNormalize, constants)
    % Non-dimensionalize or dimensionalize a state based on its type
    switch normalizationType
        case 'length'
            scale = constants.LENGTH_SCALE;
        case 'time'
            scale = constants.TIME_SCALE;
        case 'velocity'
            scale = constants.VELOCITY_SCALE;
        case 'force'
            scale = constants.VEHICLE_WEIGHT;
    end

    if flagNormalize
        % The state nees to be non-dimensionalized
        stateOut = state / scale;
    else
        % Return the state to its dimensional form
        stateOut = state * scale;
    end

