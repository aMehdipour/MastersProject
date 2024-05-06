function [CL, CD, CS] = interpAeroTables(parameters, alpha, beta)
    % Convert alpha and beta to degrees
    alphaDeg = rad2deg(alpha);
    betaDeg = rad2deg(beta);
    
    % Define the range of AoA and sideslip angles in the tables
    alphaRange = linspace(-20, 20, 5);
    betaRange = linspace(-20, 20, 5);
    
    % Calculate the polynomial fit coefficients for each alpha value
    XL = NaN(length(alphaRange), 1);
    YL = NaN(length(alphaRange), 1);
    ZL = NaN(length(alphaRange), 1);
    
    XD = NaN(length(alphaRange), 1);
    YD = NaN(length(alphaRange), 1);
    ZD = NaN(length(alphaRange), 1);
    
    YS = NaN(length(alphaRange), 1);
    ZS = NaN(length(alphaRange), 1);
    
    for i = 1:length(alphaRange)
        CLrow = parameters.CLtable(i, :);
        CDrow = parameters.CDtable(i, :);
        CSrow = parameters.CStable(i, :);
        
        % Fit second-order polynomial to CL and CD data
        polyCoeffsL = polyfit(betaRange, CLrow, 2);
        polyCoeffsD = polyfit(betaRange, CDrow, 2);
        
        XL(i) = polyCoeffsL(1);
        YL(i) = polyCoeffsL(2);
        ZL(i) = polyCoeffsL(3);
        
        XD(i) = polyCoeffsD(1);
        YD(i) = polyCoeffsD(2);
        ZD(i) = polyCoeffsD(3);
        
        % Fit first-order polynomial to CS data
        polyCoeffsS = polyfit(betaRange, CSrow, 1);
        
        YS(i) = polyCoeffsS(1);
        ZS(i) = polyCoeffsS(2);
    end
    
    % Interpolate polynomial coefficients based on alpha
    XLalpha = interp1(alphaRange, XL, alphaDeg);
    YLalpha = interp1(alphaRange, YL, alphaDeg);
    ZLalpha = interp1(alphaRange, ZL, alphaDeg);
    
    XDalpha = interp1(alphaRange, XD, alphaDeg);
    YDalpha = interp1(alphaRange, YD, alphaDeg);
    ZDalpha = interp1(alphaRange, ZD, alphaDeg);
    
    YSalpha = interp1(alphaRange, YS, alphaDeg);
    ZSalpha = interp1(alphaRange, ZS, alphaDeg);
    
    % Calculate CL, CD, and CS using equations 15-17
    CL = XLalpha * betaDeg^2 + YLalpha * betaDeg + ZLalpha;
    CD = XDalpha * betaDeg^2 + YDalpha * betaDeg + ZDalpha;
    CS = YSalpha * betaDeg + ZSalpha;
end
