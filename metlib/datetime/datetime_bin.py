#!/usr/bin/env python
# -*- coding:utf-8 -*-

# datetime_bin.py
"""This module provides functions on binning and averaging data"""

import os, sys
#import re
from datetime import datetime, timedelta
import numpy as np
#import scipy as sp
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
from .parser import *

__all__ = ['datetime_bin', 'datetime_period', 'datetime_group', 'datetime_neighbor'] 
def datetime_bin(datetimes, tdelta, starttime=None, endtime=None, return_bin_info=False):
    """This function partition a serie of datetime objects into equal timedelta bin. 
    Returns a list of np.where style tuples. If return_bin_info is True, also returns a list of (bin_start, bin_end) tuples.
    Parameters:
    datetimes is a seq of datetime objects, 
    tdelta is a timedelta object,
    starttime is a datetime object or None. When it's None, use the first datetime in the input sequence as starttime
    endtime is a datetime object or None. When it's None, use the last datetime in the input sequence as endtime
    """
    if len(datetimes) == 0:
        if return_bin_info:
            return [], []
        else:
            return []
    datetimes = parse_datetime(datetimes)
    if starttime is None:
        starttime = datetimes[0]
    else:
        starttime = parse_datetime(starttime)
    if endtime is None:
        endtime = datetimes[-1] # + tdelta / 2
    else:
        endtime = parse_datetime(endtime)
    if tdelta is None:
        tdelta = endtime - starttime
    else:
        tdelta = parse_timedelta(tdelta)

    tmpbin_start = starttime
    result = []
    bin_info = []
    tmpbin = []
    tmpbin_end = tmpbin_start + tdelta
    i = 0
    while tmpbin_start < endtime:
        if i == len(datetimes) or datetimes[i] >= tmpbin_end:
            bin_info.append((tmpbin_start, tmpbin_end))
            result.append((np.array(tmpbin, dtype='i8'), ))
            tmpbin = []
            tmpbin_start, tmpbin_end = tmpbin_end, tmpbin_end + tdelta
        elif datetimes[i] < tmpbin_start:
            i += 1
        else:
            tmpbin.append(i)
            i += 1
            
    if return_bin_info is True:
        return result, bin_info
    else:
        return result

def datetime_period(datetimes, splitpoints, starttime=None, endtime=None, return_bin_info=False, include_period_before_first_splitpoint=False):
    """This function partition a serie of datetime objects into periods according to splitpoints. 
    Returns a list of np.where style tuples. If return_bin_info is True, also returns a list of (bin_start, bin_end) tuples.
    Parameters:
    datetimes is a seq of datetime objects, 
    splitpoints is a seq of datetime objects as splitpoints,
    starttime is a datetime object or None. When it's None, use the first datetime in the input sequence as starttime
    endtime is a datetime object or None. When it's None, use the last datetime in the input sequence as endtime
    return_bin_info: default False
    include_period_before_first_splitpoint: default False
    """
    if len(datetimes) == 0:
        if return_bin_info:
            return [], []
        else:
            return []
    datetimes = parse_datetime(datetimes)
    min_dts = np.min(datetimes)
    max_dts = np.max(datetimes)
    if starttime is None:
        starttime = min_dts
    else:
        starttime = parse_datetime(starttime)
    if endtime is None:
        endtime = max_dts
    else:
        endtime = parse_datetime(endtime)
    if starttime > max_dts or endtime < min_dts:
        if return_bin_info:
            return [], []
        else:
            return []
    
    if not include_period_before_first_splitpoint:
        begs = list(splitpoints)
        ends = list(splitpoints)[1:] + [datetime.max]
    else:
        begs = [datetime.min] + list(splitpoints)
        ends = list(splitpoints) + [datetime.max]
    
#    start_i = np.argmax(splitpoints >= starttime)
#    end_i = np.argmin(splitpoints >= endtime)
# TODO: add actual starttime and endtime support
    bes = zip(begs, ends)
    periods = []
    for b, e in bes:
        periods.append( np.where((datetimes >= b) & (datetimes < e)) )

    if return_bin_info:
        return periods, bes
    else:
        return periods


def datetime_group(datetimes, tdelta_threshold, starttime=None, endtime=None, return_bin_info=False):
    """This function partition a serie of datetime objects into groups (clusters) according to a tdelta_threshold. 
    Returns a list of np.where style tuples. If return_bin_info is True, also returns a list of (bin_start, bin_end) tuples.
    Parameters:
    datetimes is a seq of datetime objects, which should be sorted.
    tdelta_threshold: a timedelta object.
    starttime is a datetime object or None. When it's None, use the first datetime in the input sequence as starttime. Not Implemented
    endtime is a datetime object or None. When it's None, use the last datetime in the input sequence as endtime. NOt Implemented
    return_bin_info: default False
    """
    if len(datetimes) == 0:
        if return_bin_info:
            return [], []
        else:
            return []
    datetimes = parse_datetime(datetimes)
    min_dts = np.min(datetimes)
    max_dts = np.max(datetimes)
    if starttime is None:
        starttime = min_dts
    else:
        starttime = parse_datetime(starttime)
    if endtime is None:
        endtime = max_dts
    else:
        endtime = parse_datetime(endtime)
    if starttime > max_dts or endtime < min_dts:
        if return_bin_info:
            return [], []
        else:
            return []
    
    basket = [0]
    pool = []
    for i in range (1, len(datetimes)):
        if datetimes[i] - datetimes[i-1] < tdelta_threshole:
            basket.append(i)
        else:
            pool.append(basket)
            basket = [i]
    pool.append(basket)
    result = [(b, ) for b in pool]
    if return_bin_info:
        return result, [(datetimes[b[0]], datetimes[b[-1]]) for b in pool] 
    else:
        return result

def datetime_neighbor(datetimes, target_datetimes, tdelta_threshold, starttime=None, endtime=None, return_bin_info=False):
    """This function finds proper datetimes which are near each target datetime point, within tdelta_threshold, from a serie of datetimes . 
    Returns a list of np.where style tuples. If return_bin_info is True, also returns a list of (bin_start, bin_end) tuples.
    Parameters:
    datetimes is a seq of datetime objects, to be choosen from.
    target_datetimes: for each one in this, find its neighborhood.
    tdelta_threshold: a timedelta object.
    starttime is a datetime object or None. When it's None, use the first datetime in the input sequence as starttime. Not Implemented
    endtime is a datetime object or None. When it's None, use the last datetime in the input sequence as endtime. NOt Implemented
    return_bin_info: default False
    """
    if len(datetimes) == 0:
        if return_bin_info:
            return [], []
        else:
            return []
    datetimes = parse_datetime(datetimes)
    min_dts = np.min(datetimes)
    max_dts = np.max(datetimes)
    if starttime is None:
        starttime = min_dts
    else:
        starttime = parse_datetime(starttime)
    if endtime is None:
        endtime = max_dts
    else:
        endtime = parse_datetime(endtime)
    if starttime > max_dts or endtime < min_dts:
        if return_bin_info:
            return [], []
        else:
            return []
    
    result = []
    bes = []
    for target in target_datetimes:
        b = target - tdelta_threshold
        e = target + tdelta_threshold
        w = np.where((datetimes>=b) & (datetimes<=e))
        result.append(w)
        bes.append((b,e))

    if return_bin_info:
        return result, bes
    else:
        return result

if __name__ == '__main__':
    pass
