import numpy as np
EARTH_RADIUS = 6378000 # m
G = np.pi**2

def atmosphere_model(altitude: float) -> float:
    """
    Atmospheric model based on https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html.
    Calculated the pressure, temperature, and density of the air at a given altitude.

    Args:
        altitude (float): Altitude above sea level in meters

    Returns:
        temp (float): Temperature in degrees Kelvin
        pressure (float): Pressure in kPa
        rhos (float): Density in kg/m^3
    """
    if altitude >= 25000:
        temp = -131.21 + 0.00299 * altitude
        pressure = 2.488 * ((temp + 273.1) / 216.6)**(-11.388)
        rho = pressure / (.2869 * (temp + 273.1))
    elif altitude < 25000 and altitude > 11000:
        temp = -56.46
        pressure = 22.65 * np.exp(1.73 - 0.000157 * altitude)
        rho = pressure / (.2869 * (temp + 273.1))
    elif altitude <= 11000:
        temp = 15.04 - 0.00649 * altitude
        pressure = 101.29 *((temp + 273.1) / 288.08)**5.256
        rho = pressure / (.2869 * (temp + 273.1))
    return temp + 273.1, pressure, rho

def get_gravity(altitude: float, earthRadius=EARTH_RADIUS, g0=G) -> float:
    g = g0 * (earthRadius / (earthRadius + altitude))
    return g
