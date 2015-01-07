#!/usr/bin/env python

import numpy as np
from .constants import *
from .humidity import *

def T2Tv(T, hum, humidity_type='RH', p=p_std, water_surface=True):
    "T to Tv. T: deg C"
    if humidity_type == 'e':
        e = hum
    elif humidity_type == 'RH':
        e = hum_rh2e(hum, T, water_surface)
    elif humidity_type == 'q':
        e = hum_q2e(hum, p)
    else:
        pass
    return (T+273.15) / (1.0-e/p*(1.0-epsilon_vapor)) - 273.15

def T2theta(T, p, p0=1000.0):
    "T to theta. T: deg C"
    return (T+273.15) * np.power(p0/p, 0.286) - 273.15
