#!/usr/bin/env python

# filter_cloud.py

#import os, sys
#import re
#from datetime import datetime, timedelta
import numpy as np
#import scipy as sp
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
from .lidar import LidarDataset

__all__ = ['filter_cloud']

def filter_cloud(data, gate, repeat=1, start_index=0, fillvalue=np.nan):
    """
    Returns: filtered_data, cloud_index
    """
    if type(data) is LidarDataset:
        data = data['data']
    gate = np.asarray(gate)
    if gate.shape != tuple() and gate.shape[-1] != 1:
        gate = gate[...,np.newaxis]
    
    cr = data > gate
    cr[..., :start_index] = False
    det = cr[..., repeat-1:].copy()
    for i in range(1,repeat):
        det = det & cr[..., repeat-1-i:-i]
    # treat as cloud if any channel has cloud
    final = np.logical_or.accumulate(det, axis=1)
    if start_index < repeat :
        start_index = repeat 
    # The trick of including one pixel under start_index makes that if there's no True value, 
    # it will return the index below start_index, which is filtered by the tricky compare below. Using start_index-1 to cope with the situation when all False(no cloud), then in the following check: if index[i] >= start_index..., start_index-1 will be filter out. Too tricky anyway.
    indexs = (final[..., start_index-1:].argmax(axis=-1) + start_index - 1).reshape(final.shape[0])
#    print indexs
    
    newdata = data.copy()
    for i in range(data.shape[0]):
        if indexs[i] >= start_index:
            newdata[i, :, indexs[i]:] = fillvalue
        else:
            indexs[i] = -1
    return newdata, indexs

if __name__ == '__main__':
    pass
