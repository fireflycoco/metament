#!/usr/bin/env python

import numpy as np

__all__ = ['exp_interp']

def exp_interp(x, x1, x2, y1, y2):
    """exp_interp interpolates in form of y = bx^{-a} , as used in interpolating AOD(Aerosol Optical Depth) from 2 wavelength.
    Parameters:
        x: the target AOD wavelength.
        x1, x2: the 2 available AOD wavelength.
        y1, y2: AOD in the 2 wavelength.
        x, x1, x2 should be scaler. y1, y2 should be in the same shape. 
    Returns:
        interpolated y in the same shape as y1 or y2.
"""
    malpha = np.log(y2 / y1) / np.log(x2 / x1)
    beta = np.exp(np.log(y1) - malpha * np.log(x1))
    y = beta * np.power(x, malpha)
    return y

