#!/usr/bin/env python2.7

# util.py

import os, sys
import re
#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
#import pandas as pd
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.ticker import FuncFormatter, FixedFormatter, ScalarFormatter, LogFormatter, FormatStrFormatter
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset
import itertools
import warnings
from metlib.misc.datatype import Null, isseq

__all__ = ['format_ticks', 'flip_legend']

def format_ticks(ax=None, fmt=None, axis='x', **kwargs):
    """Format ticks.
    ax: ax, default gca().
    fmt: date format like "%Y%m%d", or format string like "%.2f", or a func, or sequence of strs, or sequence of (position, str), default None.
    axis: 'x' | 'y' | 'z', default 'x'.
    kwargs: anything that can be applied to ticklabels, e.g., fontsize, fontweight, etc.
    """
    if not ax:
        ax = plt.gca()

    if axis == 'x':
        theaxis = ax.xaxis
    elif axis == 'y':
        theaxis = ax.yaxis
    elif theaxis == 'z':
        theaxis = ax.zaxis
    else:
        theaxis = Null

    if fmt:
        datefmt_pattern = r'%H|%I|%M|%S|%Y|%y|%m|%w|%W|%a|%A|%b|%B|%p|%z|%Z|%j|%U|%c|%x|%X'
        str_pattern = r'%'
        if callable(fmt):
            fmt = FuncFormatter(fmt)
        elif isinstance(fmt, (str, unicode, np.string_)):
            if re.search(datefmt_pattern, fmt):
                fmter = DateFormatter(fmt)
            elif re.search(str_pattern, fmt):
                fmter = FormatStrFormatter(fmt)
            else:
                fmter = Null
        elif isseq(fmt):
            if np.ndim(fmt) == 1:
                fmter = FixedFormatter(fmt)
            elif np.ndim(fmt) == 2 and np.shape(fmt)[1] == 2:
                fmt = np.array(fmt)
                pos = fmt[:, 0]
                fmter = FixedFormatter(fmt[:,1])
                theaxis.set_ticks(pos)
            else:
                fmter = Null
        else:
            fmt = Null

        if fmter:
            theaxis.set_major_formatter(fmter)

    ticklabels = theaxis.get_ticklabels()
    for key, value in kwargs.iteritems():
        little_kwargs = {key:value}
        try:
            plt.setp(ticklabels, **little_kwargs)
        except Exception as e:
            warnings.warn(e)
            

def flip_legend(ax, *args, **kwargs):
    """flip_legend() plots the legend of an axes as plt.legend(), but the labels fills rows first.
    According to Avaris' solution on stackoverflow.com:
    http://stackoverflow.com/questions/10101141/matplotlib-legend-add-items-across-columns-instead-of-down
    """
    def flip(items, ncol):
        return itertools.chain(*[items[i::ncol] for i in range(ncol)])
    
    ncol = kwargs['ncol'] if 'ncol' in kwargs else 1
    handles, labels = ax.get_legend_handles_labels()
    flipped_handles = flip(handles, ncol)
    flipped_labels = flip(labels, ncol)
    ax.legend(flipped_handles, flipped_labels, *args, **kwargs)

if __name__ == '__main__':
    pass
