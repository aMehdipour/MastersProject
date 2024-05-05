function [trimAoA, trimCD, trimCL] = calculateTrimAero(parameters, constants, velocity, altitude)
    AoARange = linspace(20, -20, 5);

    % Select the column corresponding to 0 sideslip from CL table
    CL_at_zero_sideslip = parameters.CLtable(:,3); % Assuming the 3rd column corresponds to 0 sideslip

    % Calculate weight
    weight = constants.VEHICLE_MASS * constants.GRAVITY_NOMINAL;

    % Use interp1 to find the AoA that gives the required lift
    AoA_fine = linspace(20, -20, 1000); % Fine AoA range for better resolution
    CL_interp = interp1(AoARange, CL_at_zero_sideslip, AoA_fine, 'linear');

    % Find the AoA where lift equals weight (trim AoA)
    [rho, ~, ~] = atmosphereModel(altitude);
    Lift = 0.5 * rho * (velocity^2) * CL_interp * constants.FRONTAL_AREA;
    [~, idx] = min(abs(Lift - weight)); % Find the index of the closest value to weight
    trimAoA = AoA_fine(idx);
    trimCL = CL_interp(idx);  % Also return trimCL

    % Find the drag coefficient that corresponds to the trim AoA
    trimCD = interp1(AoARange, parameters.CDtable(:,3), trimAoA, 'linear'); % Interpolate CD at trimAoA
end

