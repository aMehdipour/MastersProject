function deadbandWidth = calculateDeadbandWidth(velocity, constants)
    % Calculate the deadband width (Eq. 23 in the paper)
    deadbandWidth = constants.C1 * velocity + constants.C0;
end
