% function [trimAoA, trimCD, trimCL] = calculateTrimAero(parameters, constants, velocity, altitude)
function [maxLD, maxCD, maxCL] = calculateTrimAero(parameters, constants, velocity, altitude)
    % AoARange = linspace(20, -20, 5);
    % 
    % % Select the column corresponding to 0 sideslip from CL table
    % CL_at_zero_sideslip = parameters.CLtable(:,3); % Assuming the 3rd column corresponds to 0 sideslip
    % 
    % % Calculate weight
    % weight = constants.VEHICLE_MASS * constants.GRAVITY_NOMINAL;
    % 
    % % Use interp1 to find the AoA that gives the required lift
    % AoA_fine = linspace(20, -20, 1000); % Fine AoA range for better resolution
    % CL_interp = interp1(AoARange, CL_at_zero_sideslip, AoA_fine, 'linear');
    % 
    % % Find the AoA where lift equals weight (trim AoA)
    % [rho, ~, ~] = atmosphereModel(altitude);
    % Lift = 0.5 * rho * (velocity^2) * CL_interp * constants.FRONTAL_AREA;
    % [~, idx] = min(abs(Lift - weight)); % Find the index of the closest value to weight
    % trimAoA = AoA_fine(idx);
    % trimCL = CL_interp(idx);  % Also return trimCL
    % 
    % % Find the drag coefficient that corresponds to the trim AoA
    % trimCD = interp1(AoARange, parameters.CDtable(:,3), trimAoA, 'linear'); % Interpolate CD at trimAoA

       % Extract the range of AoA and sideslip angles from the tables
    alphaRange = linspace(-20, 20, size(parameters.CLtable, 1));
    betaRange = linspace(-20, 20, size(parameters.CLtable, 2));
    
    % Initialize variables to store the maximum L/D and corresponding CL and CD
    maxLD = -Inf;
    maxCL = 0;
    maxCD = 0;
    
    % Iterate over all combinations of AoA and sideslip angles
    for i = 1:length(alphaRange)
        for j = 1:length(betaRange)
            CL = parameters.CLtable(i, j);
            CD = parameters.CDtable(i, j);
            
            % Calculate the L/D ratio
            LD = CL / CD;
            
            % Update the maximum L/D and corresponding CL and CD if necessary
            if LD > maxLD
                maxLD = LD;
                maxCL = CL;
                maxCD = CD;
            end
        end
    end
end

