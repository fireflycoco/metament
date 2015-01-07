#!/usr/bin/env python2.7

# earth.py

import os, sys
#import re
#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
#import pandas as pd
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset

__all__ = ['latlon_distance']

def latlon_distance(lat1, lon1, lat2, lon2, R=6371.0):
    """calculate distance between 2 lat/lon points.
    Haversine Formula.
    R: earth radius.

    The code is based on the solution given by http://stackoverflow.com/questions/15736995/how-can-i-quickly-estimate-the-distance-between-two-latitude-longitude-points.
    """
    lat1 = np.deg2rad(lat1)
    lon1 = np.deg2rad(lon1)
    lat2 = np.deg2rad(lat2)
    lon2 = np.deg2rad(lon2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    c = 2.0 * np.arcsin(np.sqrt(a))
    km = c * R
    return km


if __name__ == '__main__':
    pass
