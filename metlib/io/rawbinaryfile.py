#!/usr/bin/env python

import os
import numpy as np

class RawBinaryFile(object):
    def __init__(self, fname, dtype, col_size, row_size):
        self.fname = fname
        self.col_size = col_size
        self.row_size = row_size
        self.dtype = np.dtype(dtype)
        assert os.path.getsize(self.fname) == col_size * row_size * \
            self.dtype.itemsize

        print self.dtype

    def read(self, col_start, row_start, col_num, row_num):
        data = np.empty((row_num, col_num), dtype=self.dtype)
        skip_to_pos = row_start * self.col_size * self.dtype.itemsize
        one_row_length = self.col_size * self.dtype.itemsize
        skip_line_beginning = col_start * self.dtype.itemsize

        f = open(self.fname)
        for i in xrange(row_num):
#            print i
            f.seek(skip_to_pos + i * one_row_length + skip_line_beginning)

            tmpbytes = f.read(col_num * self.dtype.itemsize)
            tmpdatarow = np.fromstring(tmpbytes, dtype=self.dtype)
#            print tmpdatarow
            data[i] = tmpdatarow[:]
        return data
        

