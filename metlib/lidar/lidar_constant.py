#!/usr/bin/env python

# lidar_constant.py

import os, sys
#import re
#from datetime import datetime, timedelta
import numpy as np
import numpy.ma as ma
#import scipy as sp
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
from .constant import lidar_sm
from .lidarutil import height_to_index

__all__ = ['get_lidar_constant']

def get_lidar_constant(data, aod, height_range, betam, elev_angle, aod_ratio=1.0, return_detail=False):
    """get lidar constant with aod observation.
    data: a LidarDataset object.
    aod: AOD value
    height_range: 2-tuple, low and high height range of lidar data for calculating lidar constant.
    betam: betam array.
    elev_angle: lidar elevation angle in degrees.
    aod_ratio: the ratio of (aod below that height) / (total aod).
    return_detail: return LC(avered), LC(detail), heights
    """
    start_i = data['first_data_bin']
    dz = data['bin_size'] / 1000.0   # convert to km
    d = data['data'][:]
    sin_elev = np.sin(np.deg2rad(elev_angle))
    intbm = np.zeros_like(betam)
    intbm[start_i:] = np.add.accumulate(betam[start_i:]) * dz
#    print intbm
    expIntSigmam = np.exp(-2.0 * intbm * lidar_sm)
#    print expIntSigmam
    low_index, high_index = height_to_index(height_range, data, elev_angle)
#    print low_index, high_index
    CE_pool = d[:,:,low_index:high_index] / (betam[low_index:high_index] * expIntSigmam[low_index:high_index] * np.exp(-2.0 * aod * aod_ratio / sin_elev) ) 
#    print CE_pool
    C_pool = CE_pool / data['energy'][..., np.newaxis]
    lc_lines = np.array(ma.masked_invalid(C_pool).mean(axis=-1)) 
#    print lc_lines
    if return_detail:
        return lc_lines, C_pool, data['distance'][low_index:high_index] * sin_elev
    else:
        return lc_lines
    


if __name__ == '__main__':
    pass
