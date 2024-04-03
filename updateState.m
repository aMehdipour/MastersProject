% Update state variables function
function [altitude, velocity, flightPathAngle, range] = updateState(altitude, velocity, flightPathAngle, range, commandedBankAngle, constants, parameters)

    % Compute non-dimensionalized lift and drag forces
    [rho, ~, ~] = atmosphereModel(altitude);
    lift = 0.5 * rho * CL * velocity^2;
    drag = 0.5 * rho * CD * velocity^2;
    
    % Equations of motion
    velocityDot = -drag - sin(flightPathAngle);
    flightPathAngleDot = (lift * cos(commandedBankAngle) - cos(flightPathAngle)) / velocity;
    altitudeDot = velocity * sin(flightPathAngle);
    rangeDot = velocity * cos(flightPathAngle);
    
    % Update state variables
    velocity = velocity + velocityDot * dt;
    flightPathAngle = flightPathAngle + flightPathAngleDot * dt;
    altitude = altitude + altitudeDot * dt;
    range = range + rangeDot * dt;
end
