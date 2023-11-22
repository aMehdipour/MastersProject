function [L, D] = computerAeroDummy(V, rho, constants)
    CL0 = 0.5; % Base coefficient of lift
    % CD0 = 0.05; % Base coefficient of drag
    LDRatio = 4; % Lift-to-drag ratio

    L = CL0 * V^2/2 * rho; % Lift proportional to square of velocity
    D = 1.5 * V^2/2 * rho; % Drag based on lift-to-drag ratio

    S = constants.REFERENCE_LENGTH; % Reference area (square meters, example value)
    L = L * S;
    D = D * S;
end
