#!/usr/bin/env python
# anomaly.py

#import os, sys
#import re
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
#from dateutil.parser import parse
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset
from datetime_range import *
from logical import *
from misc import *
from parser import *

__all__ = ['three_month_anomaly']

def three_month_anomaly(dts, data, beg_year, beg_month, end_year, end_month, day=1, least_days=0, return_base=False):
    """three_month_anomaly calculates 3-month anomaly of given data series.
    dts: datetime seq.
    data: data seq.
    beg_year, beg_month; end_year, end_month: begin/end date.
    day: day of the beg/end date
    least_days: if the sample number in a 3-month bin is less than least_days, it's omitted.
    """
    beg_dt = datetime(beg_year, beg_month, day)
    end_dt = datetime(end_year, end_month, day)
    befbeg_dt = beg_dt - TD('1M')
    aftend_dt = end_dt + TD('2M')
    tot_yms = datetime_range(befbeg_dt, aftend_dt, '1M')
    beg_dts = tot_yms[:-3]
    mid_dts = tot_yms[1:-2]
    end_dts = tot_yms[3:]
    data_shp = np.shape(data)
    stack_base_shape = (end_year-beg_year+1, 12)
    result_base_shape = ((end_year-beg_year+1) * 12, )
    if len(data_shp) > 1:
        stack_shape = stack_base_shape + data_shp[1:]
        result_shape = result_base_shape + data_shp[1:]
    else:
        stack_shape = stack_base_shape
        result_shape = result_base_shape
    stack = np.zeros(stack_shape, dtype='f8')
    stack[:] = np.nan
    for i in range(len(beg_dts)):
        beg = beg_dts[i]
        end = end_dts[i]
        yyyy = mid_dts[i].year
        mm   = mid_dts[i].month
        sel = np.ma.masked_invalid(data[np.where(datetime_is_between(beg, end, dts))])
        stack[yyyy-beg_year, mm-1] = sel.mean(axis=0)
        if len(data_shp) > 1:
            stack[yyyy-beg_year, mm-1, sel.count(axis=0) < least_days] = np.nan
        else:
            if sel.count() < least_days:
                stack[yyyy-beg_year, mm-1] = np.nan
    basevalue = np.ma.masked_invalid(stack).mean(axis=0)
    anomaly_stack = stack - basevalue[np.newaxis, ...]
    anom_res = anomaly_stack.reshape(result_shape)[beg_month-1:beg_month-1+len(beg_dts)].data
    if return_base:
        return anom_res, basevalue
    else:
        return anom_res


