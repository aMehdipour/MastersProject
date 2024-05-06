function [CL, CD, CS] = calculateAeroCoefficients(parameters, alpha, beta)
 % Convert alpha and beta to degrees
    alphaDeg = -rad2deg(alpha);
    betaDeg = rad2deg(beta);
    
    % Define the range of AoA and sideslip angles in the tables
    alphaRange = linspace(-20, 20, 5);
    betaRange = linspace(-20, 20, 5);
    
    % Interpolate the tables to a finer resolution
    finerAlphaRange = linspace(min(alphaRange), max(alphaRange), 100);
    finerBetaRange = linspace(min(betaRange), max(betaRange), 100);
    [betaMesh, alphaMesh] = meshgrid(finerBetaRange, finerAlphaRange);
    
    CLinterp = interp2(betaRange, alphaRange, parameters.CLtable, betaMesh, alphaMesh);
    CDinterp = interp2(betaRange, alphaRange, parameters.CDtable, betaMesh, alphaMesh);
    CSinterp = interp2(betaRange, alphaRange, parameters.CStable, betaMesh, alphaMesh);
    
    % Interpolate CL, CD, and CS values based on alpha and beta
    CL = interp2(finerBetaRange, finerAlphaRange, CLinterp, betaDeg, alphaDeg);
    CD = interp2(finerBetaRange, finerAlphaRange, CDinterp, betaDeg, alphaDeg);
    CS = interp2(finerBetaRange, finerAlphaRange, CSinterp, betaDeg, alphaDeg);
end
