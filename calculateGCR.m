function [greatCircleRange] = calculateGCR(lat1, lon1, lat2, lon2, radiusEarth, mode)
    if nargin < 5
        mode = 'rad';
    end

    if contains(mode, 'deg' )
        lat1 = deg2rad(lat1);
        lat2 = deg2rad(lat2);
        lon1 = deg2rad(lon1);
        lon2 = deg2rad(lon2);
    end
    
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;

    % Use the Haversine formula
    a = sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));

    greatCircleRange = radiusEarth * c;
