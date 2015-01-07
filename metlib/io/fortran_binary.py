#!/usr/bin/env python
# -*- coding:utf-8 -*-

# fortran_binary.py

import os, sys
#import re
#from datetime import datetime, timedelta
import numpy as np
#import scipy as sp
#from scipy.io.numpyio import fread as sp_fread
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab

__all__ = ['FortranBinary']
class FortranBinary(object):
    def __init__(self, fname, dtype='f8'):
        self.dtype = np.dtype(dtype)
        if self.dtype.byteorder == '>':
            rechead_dtype = '>i4'
        else:
            rechead_dtype = 'i4'
        self.records = []
        self._f = open(fname, 'rb')
        self._f.seek(0,2)
        self._filesize = self._f.tell()
        now_pos=0
        while now_pos != self._filesize :
            recsize = np.memmap(self._f, dtype=rechead_dtype, mode='c', offset=now_pos, shape=(1,))[0]
            now_pos += 4   # rechead size
            self.records.append(np.memmap(self._f, dtype=dtype, mode='c', offset=now_pos, shape=(recsize / self.dtype.itemsize,)))
            now_pos += recsize
            now_pos += 4   # rechead size

    def __len__(self):
        return len(self.records)

    def __getitem__(self, i):
        return self.records[i]

    def __del__(self):
        self._f.close()
# TODO: add rec index, etc.
