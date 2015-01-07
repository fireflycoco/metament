#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np
from matplotlib.dates import date2num
from .parser import *

__all__ = ['datetime_interp']

def datetime_interp(dest_dts, src_dts, src_values):
    """datetime_interp interpolates data to the given datetime points(dest_dts).
    Parameters:
        dest_dts: seq of datetimes. Values are interpolates to these datetime points.
        src_dts: seq of raw data's datetimes.
        src_values: seq of raw data.
    Returns:
        array of data with the same length of dest_dts.
    """
    dest_dts = parse_datetime(dest_dts)
    src_dts = parse_datetime(src_dts)
    src_dt_nums = date2num(src_dts)
    dest_dt_nums = date2num(dest_dts)
    res = np.interp(dest_dt_nums, src_dt_nums, src_values)
    return res

