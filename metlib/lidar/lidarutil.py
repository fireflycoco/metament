#!/usr/bin/env python

# lidarutil.py

import os, sys
#import re
#from datetime import datetime, timedelta
import numpy as np
#import scipy as sp
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from .constant import lidar_sm
from metlib.lidar.lidar import LidarDataset

def height_to_index(height, data, elev_angle=90.0, use_float_index=False):
    """height_to_index convert height values to lidar data index.
    height: in m.
    data: LidarDataset | equivalent dict that contains 'distance' and 'bin_size' | (resolution, offset) tuple | distance array or list.
    """
    if isinstance(data, tuple):
        resolution = data[0] * 299792458.0 / 300000000.0
        offset = data[1]
        refdist = np.arange(8000.0) * resolution + offset
    elif isinstance(data, (dict, LidarDataset)):
        refdist = data['distance']
        bin_size = data['bin_size']
        offset = refdist[0]
    elif isinstance(data, (np.ndarray, list)):
        refdist = data
        bin_size = data[1] - data[0]
        offset = refdist[0]
    else:
        raise ValueError("Wrong 'data' type: %s\n" % type(data) + self.__doc__)
    dist = np.array(height) / np.sin(np.deg2rad(elev_angle))
    float_index = (dist - offset) / bin_size 
    if use_float_index:
        return float_index
    else:
        if np.ndim(float_index) == 0:
            if np.isfinite(float_index):
                return int(np.round(float_index))
            else:
                return -1
        else:
            float_index[np.where(~np.isfinite(float_index))] = -1.0
            return np.round(float_index).astype('i')

if __name__ == '__main__':
    pass
