#!/usr/bin/env python2.7

# one_by_one.py

import os, sys
#import re
#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
#import pandas as pd
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset
from metlib.misc.datatype import Null, isinteger

__all__ = ['OnePlotter', 'OneByOneBase', 'OneByOne', 'XOneByOne', 'YOneByOne']

class OnePlotter(object):
    def __init__(self, ax, xoffset=0, yoffset=0, xratio=1.0, yratio=1.0):
        """xratio, yratio not implemented."""
        if ax is None:
            ax = plt.gca()
        self.ax = ax
        self.xoffset = xoffset
        self.yoffset = yoffset
        self.xratio = xratio
        self.yratio = yratio

    def trans_x(self, x):
        return np.array(x) * self.xratio + self.xoffset

    def trans_y(self, y):
        return np.array(y) * self.yratio + self.yoffset

    def plot(self, x, y, *args, **kwargs):
        self.ax.plot(self.trans_x(x), self.trans_y(y), *args, **kwargs)

    def scatter(self, x, y, *args, **kwargs):
        self.ax.scatter(self.trans_x(x), self.trans_y(y), *args, **kwargs)


class OneByOneBase(object):
    def __init__(self, ax=None):
        if ax is None:
            ax = plt.gca()
        self.ax = ax
            
        self._plotter_list = []

    def __getitem__(self, key):
        if isinteger(key):
            return self._plotter_list[key]
        else:
            return Null

    def __len__(self):
        return len(self._plotter_list)

class OneByOne(OneByOneBase):
    def __init__(self, xoffsets, yoffsets, xratio=1.0, yratio=1.0, ax=None):
        OneByOneBase.__init__(self, ax=ax)
        xoffsets = np.array(xoffsets)
        yoffsets = np.array(yoffsets)
        self.xratio = xratio
        self.yratio = yratio
        for i in range(len(xoffsets)):
            self._plotter_list.append(OnePlotter(self.ax, xoffset=xoffsets[i], yoffset=yoffsets[i], xratio=self.xratio, yratio=self.yratio))

class XOneByOne(OneByOneBase):
    def __init__(self, xoffsets, xratio=1.0, ax=None):
        OneByOneBase.__init__(self, ax=ax)
        xoffsets = np.array(xoffsets)
        self.xratio = xratio
        for i in range(len(xoffsets)):
            self._plotter_list.append(OnePlotter(self.ax, xoffset=xoffsets[i], xratio=self.xratio))

class YOneByOne(OneByOneBase):
    def __init__(self, yoffsets, yratio=1.0, ax=None):
        OneByOneBase.__init__(self, ax=ax)
        yoffsets = np.array(yoffsets)
        self.yratio = yratio
        for i in range(len(yoffsets)):
            self._plotter_list.append(OnePlotter(self.ax, yoffset=yoffsets[i], yratio=self.xratio))

        
if __name__ == '__main__':
    basexs = [1.0, 3.0, 4.0]
    dataxs = np.random.rand(20)
    datays = np.arange(20)

    print dataxs, datays

    xobo = XOneByOne(basexs, xratio=0.5)
    for i in range(len(xobo)):
        print xobo[i].xoffset
        xobo[i].plot(dataxs, datays)
    
    xobo.ax.set_xlim(0, 5)
    xobo.ax.set_ylim(0, 20)
    plt.show()
