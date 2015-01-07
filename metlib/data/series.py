#!/usr/bin/env python2.7

# rolling_mean.py

import os, sys
#import re
#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
import pandas as pd
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset
from metlib.data.boundary import Limiter

__all__ = ['rolling_mean']

def rolling_mean(data, window, min_periods=None, center=True, freq=None):
    """Calculate rolling mean on data. 
    data: data serie.
    window: window size.
    min_periods: min_periods.
    center: centered rolling mean.
    freq: if data is pandas.TimeSeries, resample data to this frequency. 
    """
    window = int(window)
    if isinstance(data, pd.TimeSeries):
        is_TimeSeries = True
        datats = data
        if freq:
            datats = datats.resample(freq)
        data = np.ma.masked_invalid(datats.values)
        dts = datats.index
    else:
        is_TimeSeries = False
        data = np.ma.masked_invalid(data)
    res = np.zeros(data.shape)
    res[:] = np.nan
    lendata = np.shape(data)[0]
    LI = Limiter(0, lendata-1)
    for i in range(lendata):
        if center:
            sel = data[LI[i-window/2:i-window/2+window]]
        else:
            sel = data[LI[i-window+1:i+1]]

        if min_periods and sel.count() < min_periods:
            res[i] = np.nan
        else:
            res[i] = sel.mean()
    if is_TimeSeries:
        res = pd.Series(res, index=dts)
    return res

if __name__ == '__main__':
    pass
