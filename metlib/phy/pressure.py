#!/usr/bin/env python

import numpy as np
from .constants import *

def dz_dp(p, Tv):
    "-dz/dp. p:hPa, Tv: either C/K is OK, result: m/hPa"
    return R_d*auto2K(Tv)/(p*9.80665)

