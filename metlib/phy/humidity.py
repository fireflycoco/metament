#!/usr/bin/env python

import numpy as np
from .constants import *

def hum_r2q(r):
    "Humidity conversion r(mix ratio) to q(specific humidity)"
    return r / (1.0 + r)

def hum_q2r(q): 
    "Humidity conversion q to r(mix ratio)"
    return q / (1.0 - q)

def hum_e2r(e, p=p_std):
    "Humidity conversion e(vapor pressure) to r(mix ratio). Pressure Unit: hPa."
    return epsilon_vapor * e / (p - e)

def hum_r2e(r, p=p_std):
    "Humidity conversion r(mix ratio) to e(vapor pressure). Pressure Unit: hPa."
    return hum_q2e(hum_r2q(r), p=p)

def hum_e2q(e, p=p_std):
    "Humidity conversion e(vapor pressure) to q(specific humidity). Pressure Unit: hPa."
    return epsilon_vapor * e / (p - 0.378 * e)

def hum_q2e(q, p=p_std):
    "Humidity conversion q(specific humidity) to e(vapor pressure). Pressure Unit: hPa."
    return p * q / (epsilon_vapor + 0.378 * q)

def hum_e2rho_v(e, T=T_std):
    "Humidity conversion e(vapor pressure) to rho_v(vapor density)"
    T = auto2K(T)
    return e / (R_v * T)

def hum_rho_v2e(rho_v, T=T_std):
    "Humidity conversion rho_v (vapor density) to e (wapor pressure)"
    T = auto2K(T)
    return rho_v * R_v * T

def hum_es(T, method='magnus'):
    """Compute e_satuation for water surface. 
    method: 'magnus'
    """
    if method == 'magnus':
        return hum_es_magnus(T)
    else:
        raise ValueError("hum_es method='%s' not implemented" % method)

def hum_esi(T, method='magnus'):
    """Compute e_satuation for ice surface.
    method: 'magnus'
    """
    if method == 'magnus':
        return hum_esi_magnus(T)
    else:
        raise ValueError("hum_esi method='%s' not implemented" % method)

def hum_es_magnus(T):
    "Compute e_satuation for water surface with Magnus formula."
    T = auto2K(T)
    return 6.1078 * np.exp(17.2693882*(T-273.16)/(T-35.86))

def hum_esi_magnus(T):
    "Compute e_satuation for ice surface with Magnus formula"
    T = auto2K(T)
    return 6.1078 * np.exp(21.8745584*(T-273.16)/(T-7.66))

def hum_rh2e(rh, T, water_surface=True):
    "Humidity conversion RH to e (vapor pressure)"
    if water_surface:
        return rh / 100.0 * hum_es(T)
    else:
        return rh / 100.0 * hum_esi(T)

def hum_rh2q(rh, T, p=p_std, water_surface=True):
    "Humidity conversion RH to q (specific humidity)"
    return hum_e2q(hum_rh2e(rh, T, water_surface), p)

def hum_rh2r(rh, T, p=p_std, water_surface=True):
    "Humidity conversion RH to r (mixing ratio)"
    return hum_e2r(hum_rh2e(rh, T, water_surface), p)

def hum_e2rh(e, T, water_surface=True):
    "Humidity conversion e (vapor pressure) to RH"
    if water_surface:
        return e / hum_es(T) * 100.0
    else:
        return e / hum_esi(T) * 100.0

def hum_r2rh(r, T, p=p_std, water_surface=True):
    "Humidity conversion r (mixing ratio) to RH"
    return hum_e2rh(hum_r2e(r, p=p), T, water_surface)

def hum_q2rh(q, T, p=p_std, water_surface=True):
    "Humidity conversion q (specific humidity) to RH"
    return hum_e2rh(hum_q2e(q, p=p), T, water_surface)

def hum_tdew2rh(Tdew, T, *args, **kwargs):
    "Humidity conversion Tdew to RH"
    T = auto2C(T)
    Tdew = auto2C(Tdew)
    RH = 100.0*(np.exp((17.625*Tdew)/(243.04+Tdew))/np.exp((17.625*T)/(243.04+T)))
    return RH

def hum_rh2tdew(rh, T, *args, **kwargs):
    "Humidity conversion RH to Tdew (C)"
    T = auto2C(T)
    Tdew = 243.04*(np.log(rh/100.0)+((17.625*T)/(243.04+T)))/(17.625-np.log(rh/100.0)-((17.625*T)/(243.04+T)))
    return Tdew
