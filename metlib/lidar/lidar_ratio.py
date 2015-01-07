#!/usr/bin/env python

# get_sa.py

import os, sys
#import re
#from datetime import datetime, timedelta
import numpy as np
#import scipy as sp
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
from .fernald import fernald
from .lidarutil import height_to_index

__all__ = ['get_lidar_ratio']

def get_lidar_ratio(data, aod, lidar_constant, betam, fill_index=0, fill_aver_num=1, C_contains_E=False, search_range=(10.0, 100.0, 0.1), maxheight=5000.0, elev_angle=90.0):
    """calculate lidar_ratio from aod data.
    data: a LidarDataset object
    aod: aod
    lidar_constant: lidar constant
    betam: molecular backscatter
    fill_index, fill_aver_num: data below this index will be filled with fill_aver_num samples above
    search_range: tuple of start, end, step of lidar ratio value to be tried.

    returns an array of lidar_ratios with the same length of data.
    """
    test_sa = np.arange(*search_range)
    test_res = np.zeros( (data.dims['TIME'], data.dims['CHANNEL'], len(test_sa) ) )
    try:
        start_i = data['first_data_bin']
    except:
        start_i = 0
    max_index = height_to_index(maxheight, data, elev_angle)
    max_index = np.min((max_index, data.dims['BIN'] - 1))
    check_i_beg = height_to_index(200.0, data, elev_angle)
    check_i_end = height_to_index(300.0, data, elev_angle)
#    print check_i_beg, check_i_end
    #print max_index
    dz = data['bin_size'] * 1E-3 * np.sin(np.deg2rad(elev_angle))
    for i in range(len(test_sa)):
        sigma_a = fernald(data, lidar_constant, test_sa[i], betam, fill_index, fill_aver_num, C_contains_E=C_contains_E, apply_on_data=False)
        to_check = sigma_a[..., check_i_beg:check_i_end]
        if np.any(to_check <= 0.0) or np.any(to_check > 0.8):
#            print "Warning: %s " % ( test_sa[i] ,)
            test_res[:,:,i] = np.nan
        else:
            test_res[:,:,i] = np.nansum(sigma_a[..., start_i:max_index+1], axis=-1) * dz
    res_dist = np.abs(test_res - aod)
#    print res_dist.shape
    min_indice = np.nanargmin(res_dist, axis=-1)
#    print min_indice
#    plt.plot(res_dist[0,0,:])
#    plt.show()
    final_res = np.zeros( (data.dims['TIME'], data.dims['CHANNEL']) )
    final_res[:] = np.nan

    for t in range(data.dims['TIME']):
        for c in range(data.dims['CHANNEL']):
            the_i = min_indice[t,c]
            if np.isnan(the_i) or res_dist[t, c, the_i] >= 0.2:
                final_res[t,c] = np.nan
            else:
                final_res[t, c] = test_sa[ min_indice[t, c] ]
    return final_res

if __name__ == '__main__':
    pass
