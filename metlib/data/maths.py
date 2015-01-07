#!/usr/bin/env python2.7

# maths.py

import os, sys
#import re
#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset

__all__ = ['second_derivate', 'int_sign']

def second_derivate(sig):
    sig = np.array(sig)
    res = np.zeros_like(sig)
    res[1:-1] = sig[:-2] + sig[2:] - sig[1:-1] * 2.0
    res[0] = np.nan
    res[-1]  = np.nan
    return res

def int_sign(a):
    s = np.sign(a)
    if np.isscalar(s):
        if -np.isfinite(s):
            s = 0.0
        return int(s)
    else:
        s[np.where(-np.isfinite(s))] = 0.0
        s = s.astype('i')
        return s
