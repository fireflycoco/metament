#!/usr/bin/env python

# data_2d_bin.py
"""This module provides functions on binning and averaging data"""

import os, sys
import numpy as np

__all__ = ['data_2d_bin'] 
def data_2d_bin(data, xs, ys, x_points, y_points):
    """This function partitions a seq of data with coordinates (xs, ys), into 2d bins defined by split points (x_points, y_points). 

Parameters:
    data: a seq of data.
    xs, ys: seqs of coordinates of the data.
    x_points, y_points: seqs of split points. 
Returns:
    A 2d list with the shape of (len(y_points)-1, len(x_points)-1).
    """
    assert len(data) == len(xs) == len(ys)
    if len(data) == 0:
        return [[[]]]
    res = []
    for jy in range(len(y_points)-1):
        to_append = []
        for ix in range(len(x_points)-1):
            to_append.append([])
        res.append(to_append)
    
    for jy in range(len(y_points)-1):
        wy = np.where( (y_points[jy] <= ys) & (ys < y_points[jy+1]) )
        sub_data = data[wy]
        sub_xs = xs[wy]
        for ix in range(len(x_points)-1):
            wx = np.where( (x_points[ix] <= sub_xs) & (sub_xs < x_points[ix+1]) )
            res[jy][ix].extend(sub_data[wx])

    return res

if __name__ == '__main__':
    a = np.random.rand(100)
    xs = np.random.rand(100) * 100.0
    ys = np.random.rand(100) * 100.0
    print xs, ys
    print data_2d_bin(a, xs, ys, [0, 20, 50, 100], [24, 50, 100])

