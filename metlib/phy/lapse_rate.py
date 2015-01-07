import numpy as np

from .humidity import hum_es, hum_e2r

__all__ = ['gamma_d', 'gamma_s']

# Using np.multiply and np.add and np.square to keep scaler.
gamma_d = 9.8 # K/km

g = 9.8076 # m / s^2
Hv_0 = 2501.0E3 # J / kg
Hv_100 = 2250.0E3 # J / kg
Rsd = 287.0 # J kg^-1 K^-1
Rsw = 462.0 # J kg^-1 K^-1
cpd = 1004.07 # J kg^-1 K^-1

def gamma_s(t, p):
    """Calculate saturated adiabatic lapse rate(gamma_s).
    t: temperature in deg C.
    p: pressure in hPa.
    Returns: gamma_s in K/km
    """
    Hv = (t / 100.0) * (Hv_100 - Hv_0) + Hv_0
    es = hum_es(t)
    r = 0.622 * es / (p - es)
    k = t + 273.15
    return g * (1.0 + Hv*r/(Rsd*k)) / (cpd + np.square(Hv) * r/ (Rsw * np.square(k))) * 1000.0

