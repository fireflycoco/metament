#!/usr/bin/env python

# alias.py
import os
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

__all__ = ['ma_out']

def ma_out(arr, vmin=-np.inf, vmax=np.inf, fill_value=None):
    """masked values outside (vmin, vmax) and invalid values.
    and fill the masked part with fill_value if it is not None.
    """
    res = np.ma.masked_outside(np.ma.masked_invalid(arr), vmin, vmax)
    if fill_value is not None:
        res = res.filled(fill_value)
    return res
