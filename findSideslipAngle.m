function beta = findSideslipAngle(parameters, alpha, requiredCS)
    % Convert alpha to degrees
    alphaDeg = rad2deg(alpha);
    
    % Define the range of AoA and sideslip angles in the tables
    alphaRange = linspace(-10, 10, 5);
    betaRange = linspace(-10, 10, 5);
    
    % Interpolate the CS table to a finer resolution
    finerAlphaRange = linspace(min(alphaRange), max(alphaRange), 100);
    finerBetaRange = linspace(min(betaRange), max(betaRange), 100);
    [betaMesh, alphaMesh] = meshgrid(finerBetaRange, finerAlphaRange);
    CSinterp = interp2(betaRange, alphaRange, parameters.CStable, betaMesh, alphaMesh);
    
    % Find the row index corresponding to the closest alpha value in the finer range
    alphaIdx = findClosestIndex(finerAlphaRange, alphaDeg);
    
    % Extract the CS values corresponding to the alpha row
    CSrow = CSinterp(alphaIdx, :);
    
    % Find the sideslip angle that produces the required CS value
    beta = interp1(CSrow, finerBetaRange, requiredCS);
    
    % Convert beta to radians
    beta = deg2rad(beta);
end
