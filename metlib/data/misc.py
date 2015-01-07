#!/usr/bin/env python2.7

# rolling_mean.py

import os, sys
#import re
from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset

__all__ = ['nearest_i', 'where2slice']

def nearest_i(arr, target):
    """return the nearest index for target in arr.
    Supporting only 1D data for now.
    """
    dists = np.abs(np.array(arr) - target)
    if isinstance(dists[0], timedelta):
        dists = np.array([td.total_seconds() for td in dists])
    return np.nanargmin(dists)

def where2slice(w):
    """convert np.where() results to slice if possible.
    """
    if len(w) != 1:
        raise ValueError("where2slice supports only 1D where results, given w is %d" % len(w))

    if len(w[0]) == 0:
        return slice(0, 0, 1)
    elif len(w[0]) == 1:
        return slice(w[0][0], w[0][0] + 1, 1)
    else:
        step = w[0][1] - w[0][0]
        start = w[0][0]
        end = w[0][-1] + step
        #TODO: check match
        return slice(start, end, step)

if __name__ == '__main__':
    from metlib.kits import *
    dts = datetime_range(20140101, 20140201, '1d')
    dt = T(20140107)
    print nearest_i(dts, dt)
