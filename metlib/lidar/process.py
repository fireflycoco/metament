#!/usr/bin/env python

# process.py

import os, sys
#import re
#from datetime import datetime, timedelta
import numpy as np
#import scipy as sp
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
from lidar import LidarDataset
__all__ = ['correct_background', 'correct_afterpulse',
        'correct_overlap', 'correct_distance',
        'correct_energy',
        'fill_lower_part',
        'zero_blind_range', 
        'brutal_denoise', 'gentle_denoise', 'find_first_index_below_gate']

def correct_background(data, sample_number=None):
    """corrects background noise.
    data: a LidarDataset object
    sample_number: use the last sample_number bins as background.
        if sample_number is None, use lidar data's prvided background.
    """
    if sample_number is None:
        bg = data['background']
    elif sample_number == 0:
        return
    else:
        bg = data['data'][:,:,-sample_number:].mean(axis=-1)
    bg = bg[..., np.newaxis]
    data['data'] -= bg
    data.desc += ',background corrected'

def correct_afterpulse(data, ap_data, zero_check_max_index=30):
    """correct afterpulse noise.
    data: a LidarDataset object.
    ap_data: afterpulse data (energy normalized). 
    zero_check_max_index: set negative results below this index as 0.
    """
    min_len = np.min((data.dims['BIN'], ap_data.shape[-1]))
    zcmi = np.min((min_len, zero_check_max_index))
    d = data['data']
    d[...,:min_len] -= ap_data[..., :min_len] * data['energy'][..., np.newaxis]
    d[...,:zcmi][np.where(d[...,:zcmi] < 0.0)] = 0.0
    
    data.desc += ',afterpulse corrected'

def correct_overlap(data, ol_data):
    """correct overlap.
    data: a LidarDataset object.
    ol_data: overlap data (values > 1 in the near range). 
    """
    min_len = np.min((data.dims['BIN'], ol_data.shape[-1]))
    data['data'][...,:min_len] *= ol_data[..., :min_len]
    data.desc += ',overlap corrected'

def correct_distance(data):
    """correct distance.
    data: a LidarDataset object.
    """
    data['data'] *= (data['distance'] ** 2 * 1E-6)
    data.desc += ',distance corrected'

def correct_energy(data):
    """correct energy.
    data: a LidarDataset object.
    """
    data['data'] /= data['energy'][..., np.newaxis]
    data.var_aver_methods['data'] = 'mean'
    data.desc += ',energy corrected'

def fill_lower_part(data, index, aver_num=3):
    """fill the lower part with value above.
    data: a LidarDataset object or an array
    index: to be filled data ends here, valid data starts here.
    aver_num: use data[..., index:index+aver_num].mean as fill value.
    """
    if type(data) is LidarDataset:
        d = data['data']
    else:
        d = data
    to_fill = np.ma.masked_invalid(d[..., index:index+aver_num]).mean(axis=-1)[..., np.newaxis]
    d[..., :index] = to_fill


def zero_blind_range(data):
    """make the blind range zero.
    data: a LidarDataset object.
    """
    try:
        start_i = data['first_data_bin']
    except:
        start_i = 0
    data['data'][...,:start_i] = 0.0

def brutal_denoise(data, gate, **kwargs):
    if isinstance(data, LidarDataset):
        data = data['data']
    indice = find_first_index_below_gate(data, gate, **kwargs)
    if len(np.shape(indice)) == 2:
        for t in range(data.shape[0]):
            for c in range(data.shape[1]):
                i = indice[t,c]
                if i != -1:
                    data[t, c, i:] = gate
    elif len(np.shape(indice)) == 1:
        for t in range(data.shape[0]):
            i = indice[t]
            if i != -1:
                data[t, i:] = gate

def gentle_denoise(data, gate, gate2=None, **kwargs):
    if gate2 is None:
        gate2 = gate
    if isinstance(data, LidarDataset):
        data = data['data']
    indice = find_first_index_below_gate(data, gate, **kwargs)
    if len(np.shape(indice)) == 2:
        for t in range(data.shape[0]):
            for c in range(data.shape[1]):
                i = indice[t,c]
                if i != -1:
                    aver = np.ma.masked_invalid(data[t, c, i:]).mean()
                    if aver < gate2:
                        aver = gate2
                    data[t,c,i:] = aver
    elif len(np.shape(indice)) == 1:
        for t in range(data.shape[0]):
            i = indice[t]
            if i != -1:
                aver = np.ma.masked_invalid(data[t, i:]).mean()
                if aver < gate2:
                    aver = gate2
                data[t,i:] = aver

def find_first_index_below_gate(data, gate, start_i=0):
    """find the index of first small value in lidar data.
    data can be LidarDataset, array[TIME, CHANNEL, BIN] or array[TIME, BIN].
    gate is the threshold.
    start_i is the first value usable.
    """
    if isinstance(data, LidarDataset):
        data = data['data']
    # to treat nan as values below gate
    cre = -(data >= gate)
    cre[...,:start_i] = False
    res = cre.argmax(axis=-1)
    # for the case all cre value are false
    if start_i > 0:
        res[np.where(res==0)] = -1
    else:
        res[np.where((data[..., 0] >= gate) & (res == 0))] = -1
    return res

if __name__ == '__main__':
    pass
