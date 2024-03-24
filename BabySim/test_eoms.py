import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
g0 = 9.81  # gravitational acceleration at the Earth's surface, m/s^2
R0 = 6371000  # Earth's radius, m

def compute_aero_dummy(V, FRONTAL_AREA=0.3849):
    """Compute aerodynamic forces using a dummy model."""
    CL0 = 0.5  # Base coefficient of lift
    LDRatio = 4  # Lift-to-drag ratio
    L = CL0 * V**2 * FRONTAL_AREA  # Lift
    D = L / LDRatio  # Drag
    return L, D

def energy_relative_eoms(e, state, bank_angle):
    """Defines the energy-relative equations of motion to be integrated."""
    s, r, gamma = state  # Unpack the state vector
    V = np.sqrt(2 * (1/r - e))  # Compute velocity from specific energy and radial distance
    L, D = compute_aero_dummy(V)  # Compute aerodynamic forces
    dsde = -np.cos(gamma) / (r * D)
    drde = np.sin(gamma) / D
    dgamma_de = (1 / (D * V**2)) * (L * np.cos(bank_angle) + (V**2 - (1/r)) * np.cos(gamma))
    return [dsde, drde, dgamma_de]

# Initial conditions
V_init = 8142  # m/s
h_init = 100000  # m
gamma_init_deg = -5.5  # degrees
sigma_init_deg = 0  # degrees

# Final conditions
V_final = 690  # m/s
h_final = 25000  # m

# Nondimensionalize initial conditions
V_nd = V_init / np.sqrt(g0 * R0)
r_nd = (R0 + h_init) / R0
gamma_init_rad = np.radians(gamma_init_deg)
sigma_init_rad = np.radians(sigma_init_deg)

# Nondimensionalize final conditions
V_final_nd = V_final / np.sqrt(g0 * R0)
r_final_nd = (R0 + h_final) / R0

# Compute specific energies
e_0 = 1/r_nd - V_nd**2 / 2
e_f = 1/r_final_nd - V_final_nd**2 / 2

# Initial state vector for integration [s, r, gamma]
initial_state = [1000, r_nd, gamma_init_rad]  # Assuming s starts from 0

# Integrate the EOMs
solution = solve_ivp(fun=lambda e, y: energy_relative_eoms(e, y, sigma_init_rad),
                     t_span=(e_0, e_f), 
                     y0=initial_state, 
                     method='RK45')

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(solution.t, solution.y[0], label='Range-to-go')
plt.plot(solution.t, solution.y[1], label='Radial distance')
plt.plot(solution.t, solution.y[2], label='Flight path angle')
plt.legend()
plt.xlabel('Specific Energy (nondimensional)')
plt.ylabel('State Variables')
plt.title('Integration of Energy-Relative EOMs')
plt.grid(True)
plt.show()

