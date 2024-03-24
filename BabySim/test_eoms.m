% Constants
g0 = 9.81;  % gravitational acceleration at the Earth's surface, m/s^2
R0 = 6371000;  % Earth's radius, m
FRONTAL_AREA = 0.3849;  % Default frontal area

% Initial conditions
V_init = 8142;  % m/s
h_init = 100000;  % m
gamma_init_deg = -5.5;  % degrees
sigma_init_deg = 0;  % degrees

% Final conditions
V_final = 690;  % m/s
h_final = 25000;  % m

% Nondimensionalize initial conditions
V_nd = V_init / sqrt(g0 * R0);
r_nd = (R0 + h_init) / R0;
gamma_init_rad = deg2rad(gamma_init_deg);
sigma_init_rad = deg2rad(sigma_init_deg);

% Nondimensionalize final conditions
V_final_nd = V_final / sqrt(g0 * R0);
r_final_nd = (R0 + h_final) / R0;

% Compute specific energies
global e_0 e_f
e_0 = 1/r_nd - V_nd^2 / 2;
e_f = 1/r_final_nd - V_final_nd^2 / 2;

% Options to stop integration when the result is complex
options = odeset('Events', @eventsFcn);

% Integrate the EOMs
[E, Y] = ode45(@(e, y) energy_relative_eoms(e, y, sigma_init_rad), [e_0 e_f], [1000, r_nd, gamma_init_rad], options);

% Storing diagnostic data (lift, drag, velocity)
L_values = zeros(length(E), 1);
D_values = zeros(length(E), 1);
V_values = zeros(length(E), 1);
for i = 1:length(E)
    e = E(i);
    r = Y(i, 2);
    V = sqrt(2 * (1/r - e)) * sqrt(g0 * R0);  % Convert back to dimensional
    [L, D] = compute_aero_dummy(V);
    L_values(i) = L;
    D_values(i) = D;
    V_values(i) = V;
end

% Convert the nondimensional values back to dimensional for plotting
S_dimensional = Y(:,1) * R0;  % Convert range-to-go back to meters
R_dimensional = (Y(:,2) - 1) * R0;  % Convert radial distance back to altitude above the surface
Gamma_dimensional = rad2deg(Y(:,3));  % Convert flight path angle back to degrees

% Plotting the results on separate charts
figure;
subplot(4,2,1);
plot(E, S_dimensional, 'DisplayName', 'Range-to-go');
xlabel('Specific Energy (nondimensional)');
ylabel('Range-to-go (m)');
title('Range-to-go vs Specific Energy');
grid on;

subplot(4,2,2);
plot(E, R_dimensional, 'DisplayName', 'Altitude');
xlabel('Specific Energy (nondimensional)');
ylabel('Altitude (m)');
title('Altitude vs Specific Energy');
grid on;

subplot(4,2,3);
plot(E, Gamma_dimensional, 'DisplayName', 'Flight Path Angle');
xlabel('Specific Energy (nondimensional)');
ylabel('Flight Path Angle (degrees)');
title('Flight Path Angle vs Specific Energy');
grid on;

subplot(4,2,4);
plot(E, L_values, 'DisplayName', 'Lift');
xlabel('Specific Energy (nondimensional)');
ylabel('Lift (N)');
title('Lift vs Specific Energy');
grid on;

subplot(4,2,5);
plot(E, D_values, 'DisplayName', 'Drag');
xlabel('Specific Energy (nondimensional)');
ylabel('Drag (N)');
title('Drag vs Specific Energy');
grid on;

subplot(4,2,6);
plot(E, V_values, 'DisplayName', 'Velocity');
xlabel('Specific Energy (nondimensional)');
ylabel('Velocity (m/s)');
title('Velocity vs Specific Energy');
grid on;

% Marking the point where the solution might become imaginary
if ~isempty(E)
    % Use the last computed point before stopping as the potential imaginary threshold
    xline(E(end), 'r--', 'Label', 'Potential Imaginary Threshold', 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'middle');
end

% Function definitions
function dyde = energy_relative_eoms(e, y, bank_angle)
    s = y(1);
    r = y(2);
    gamma = y(3);
    
    sqrt_argument = 2 * (1/r - e);

    V = sqrt(sqrt_argument);
    
    [L, D] = compute_aero_dummy(V);
    dsde = -cos(gamma) / (r * D);
    drde = sin(gamma) / D;
    dgamma_de = 0; % (1 / (D * V^2)) * (L * cos(bank_angle) + (V^2 - (1/r)) * cos(gamma));
    dyde = [dsde; drde; dgamma_de];
end

function [L, D] = compute_aero_dummy(V)
    FRONTAL_AREA = 0.3849;  % Define FRONTAL_AREA within the function
    CL0 = -0.01;  % Base coefficient of lift
    CD0 = 1.5;
    LDRatio = 0.5;  % Lift-to-drag ratio
    L = CL0 * V^2 * FRONTAL_AREA;  % Lift
    D = CD0 * V^2 * FRONTAL_AREA;  % Drag
end


function [value, isterminal, direction] = eventsFcn(~, y)
    global e_0
    r = y(2);
    sqrt_argument = 2 * (1/r - e_0);  % Use e_0 as reference
    
    value = double(sqrt_argument >= 0);  % Stop if sqrt_argument < 0
    isterminal = 1;  % Stop the integration
    direction = 0;  % All zeros are to be located
end
