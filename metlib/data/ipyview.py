#!/usr/bin/env python2.7

# ipyview.py

import os, sys
#import re
#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
from matplotlib import mlab
import pandas as pd
from netCDF4 import Dataset
from metlib.kits import *
from metlib.color import *
from metlib.data.series import *
from metlib.data.misc import *
from metlib.data.boundary import *
from metlib.lidar import *

__all__ = ['ipyview']

speaker = sys.stdout
warner = sys.stderr
globald = sys.modules['__main__'].__dict__

def load_everything(fname, basket_dest=''):
    try:
        fname = expand_path(fname)
        ext = get_ext(fname)
        if ext == '.npy':
            return np.load(fname), "npy"
        elif ext in ('.nc', '.ncf'):
            data = Dataset(fname, mode='r')
            if 'lidarname' in data.ncattrs():
                return LidarDataset(fname), "LidarDataset"
            else:
                return data, "netcdf"
        elif ext in ('.h5', '.h5f', '.hdf', '.hdf5'):
            return pd.HDFStore(fname, mode='r'), "pd.HDFStore"
        elif ext in ('.csv'):
            return pd.DataFrame.from_csv(fname), "pd.DataFrame"
        elif ext in ('.zip'):
            if basket_dest:
                globald[basket_dest] = dict()
                varnames = loadbasket(fname, dest=globald[basket_dest])
            else:
                varnames = loadbasket(fname)
            return varnames, "basket"
        elif ext in ('.pickle', '.pic'):
            return loadpickle(fname), "pickle"
        elif ext in ('.txt'):
            return np.loadtxt(fname), "txt"
    except Exception as e:
        warner.write("Error while loading : %s \n" % fname)
        warner.write(e)
        warner.write('\n')

def ipyview(*args, **kwargs):
    """ipyview automatically loads data files as variables in global scope.
    args: filenames. Loaded varnames will be data or data1, data2 ...
    kwargs: filenames for variables with specified names, e.g., a='a.npy', b='some.csv'.
"""

    if len(args) == 1:
        loaded_data, loaded_type = load_everything(args[0])
        print "----------------"
        if loaded_type == "basket":
            print "Loaded Basket variables:", loaded_data
        else:
            data = loaded_data
            globald['data'] = data
            print "file: %s, recognized as: %s" % (args[0], loaded_type)
            print ">>>", "data :"
            print data
    else:
        for i, fname in enumerate(args):
            print "----------------"
            varname = 'data%d' % i
            loaded_data, loaded_type = load_everything(fname)
            if loaded_type == "basket":
                print "Loaded Basket variables:", loaded_data
            else:
                globald[varname] = loaded_data
                print "file: %s, recognized as: %s" % (fname, loaded_type)
                print ">>>", varname, ":"
                print globald[varname]

    for varname, fname in kwargs.iteritems():
        print "----------------"
        loaded_data, loaded_type = load_everything(fname, basket_dest=varname)
        if loaded_type == "basket":
            print "Loaded Basket file", fname, "into dict:", varname, ", which contains:", loaded_data
        else:
            globald[varname] = loaded_data
            print "file: %s, recognized as: %s" % (fname, loaded_type)
            print ">>>", varname, ":"
            print globald[varname]

if __name__ == '__main__':
    pass
