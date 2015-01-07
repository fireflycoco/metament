#!/usr/bin/env python2.7

# peak.py

import os, sys
import re
#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset
from metlib.misc.maths import int_sign, second_derivate

__all__ = ['Peak', 'locate_peak', 'measure_peak', 'peak_pattern', 'parse_peak', 'parse_peaks', 'load_peaks_list', 'save_peaks_list' , 'peak2str']

class Peak(object):
    def __init__(self, center, lower=None, upper=None, depth=0.0, volume=0.0):
        self.center = center
        self.lower = lower if lower is not None else center
        self.upper = upper if upper is not None else center
        self.depth = depth
        self.volume = volume
        self.sign = int_sign(self.depth)
        
    def __repr__(self):
        return "Peak<%d,%d,%d,%E,%E>" % (self.center, self.lower, self.upper, self.depth, self.volume)

peak_pattern = r'Peak<([^<>\s]*)[,;\s]+([^<>\s]*)[,;\s]+([^<>\s]*)[,;\s]([^<>\s]*)[,;\s]+([^<>\s]*)>|Peak<>'
peak_re = re.compile(peak_pattern)

def _match2peak(mstrs):
    if len(mstrs) == 5:
        if mstrs[0] == '' or mstrs[0] == None:
            p = None
        else:
            p = Peak(center=int(mstrs[0]), 
                    lower=int(mstrs[1]), 
                    upper=int(mstrs[2]),
                    depth=float(mstrs[3]),
                    volume=float(mstrs[4]) )
    else:
        p = None
    return p

def peak2str(peak):
    if peak is None:
        return 'Peak<>'
    else:
        return repr(peak)

def parse_peak(s):
    m = re.search(peak_re, s)
    if m:
        mstrs = m.groups()
        return _match2peak(mstrs)
    else:
        return None

def parse_peaks(s):
    mstrs_list = re.findall(peak_re, s)
    res = []
    for mstrs in mstrs_list:
        res.append(_match2peak(mstrs))
    return res

def save_peaks_list(fname, peaks_list):
    f = open(fname, 'w')
    for peaks in peaks_list:
        for p in peaks:
            f.write(peak2str(p))
            f.write(' ')
        f.write('\n')
    f.close()

def load_peaks_list(fname):
    peaks_list = []
    f = open(fname)
    for line in f:
        if line.lstrip().startswith('#'):
            continue
        peaks = parse_peaks(line)
        peaks_list.append(peaks)
    f.close()
    return peaks_list

def locate_peak(sig):
    """
    Return
    ------
    an array with the same shape as sig, where 0: not peak or valley; 1/-1: peak/valley; 2/-2: pos valley/neg peak.
    """
    sig = np.array(sig)
    cr1 = np.zeros_like(sig, dtype='i')
    cr2 = np.zeros_like(sig, dtype='i')
    cr1[1:] = int_sign(sig[1:] - sig[:-1])
    cr2[:-1]  = int_sign(sig[:-1] - sig[1:])
    sigsign = int_sign(sig)
    cr = cr1 * cr2
    cr[cr <= 0] = 0
    cr[np.where((cr==1) & (cr1 * sigsign == -1))] = 2
    cr[np.where((cr==1) & (cr2 * sigsign == -1))] = 2
#    cr[cr==1] = 2
#    print cr
#    print sigsign
    cr *= sigsign
    return cr

def measure_peak(sig, use_inflection=True, return_allinfo=False):
    """measure_peak measures peaks from 1D/2D signal (2D signal are regarded as columns of individual signal lines), with all the information for the Peak class, e.g., peak center position (indices), peak size (lower/upper position), peak depth and peak volume.
    
    Parameters
    ----------
    sig:
    use_inflection: whether to use inflection points to split peaks.
    return_allinfo: whether to return raw peak array, crossing zero point array and inflection point array. 

    Return
    ------
    list/lists of Peak. And if return_allinfo is True, also returns 3 more arrays: raw peaks, crossing zero points and inflection points.
    
    For all the arrays:
        0: not anything.
    For raw peak array:
        1/-1: peak/valley; 2/-2: pos valley/neg peak; 
    For crossing zero array and inflection point array:
        4/-4: pos/neg crossing zero; 8/-8: pos/neg inflection point.
    """
    sig = np.array(sig)
    cr = locate_peak(sig)
    cr_crosszero = np.zeros_like(cr)
    cr_inflection = np.zeros_like(cr)

    # cross zero points
    cr_cr1 = -int_sign(sig[1:] * sig[:-1])
    cr_cr2 = -int_sign(sig[:-1] * sig[1:])
    cr_cr1[cr_cr1<0] = 0
    cr_cr2[cr_cr2<0] = 0
    cr_crosszero[1:] = cr_cr1
    cr_crosszero[:-1] += cr_cr2
    cr_crosszero = int_sign(cr_crosszero * sig) * 4

    # inflection points
    d2 = second_derivate(sig)
    d2p = locate_peak(d2)
    d2p[np.where( np.abs(d2p) != 1 )] = 0
    d2p[np.where( ((d2p==1) & (sig<0)) | ((d2p==-1) & (sig>0)) )] = 0
    cr_inflection[np.where(d2p==-1)] = 8
    cr_inflection[np.where(d2p==1)] = -8
   
    if use_inflection:
        cr_combine = cr + cr_inflection + cr_crosszero 
    else:
        cr_combine = cr + cr_crosszero

    oned = False
    if len(np.shape(sig)) == 1:
        oned = True
        sig = sig[:, np.newaxis]
   
    peaks_list = []
    for i in range(np.shape(sig)[1]):
        pvs = np.where(np.abs(cr[:,i]) == 1)[0]
        lims = np.where(np.abs(cr_combine[:,i]) >= 2)[0]
        if len(pvs) == 0 :
            peaks_list.append([])
            continue
        if np.shape(lims)[0] == 0:
            lower_pos = pvs
            upper_pos = pvs
        else:
            lower_arr = (pvs > lims[:, np.newaxis])
            upper_arr = (pvs < lims[:, np.newaxis])
            lower_arr_r = np.flipud(lower_arr)
            upper_pos_i = np.argmax(upper_arr, axis=0)
            upper_pos = lims[(upper_pos_i, )]
            w_upper_none = np.where(upper_arr[-1,:] == False)
            upper_pos[w_upper_none] = pvs[w_upper_none]
            lower_pos_r_i = np.argmax(lower_arr_r, axis=0)
            lower_pos_i = len(lims) - 1 - lower_pos_r_i
            lower_pos = lims[(lower_pos_i, )]
            w_lower_none = np.where(lower_arr[0, :] == False)
            lower_pos[w_lower_none] = 0

        peaks = []
        for center, lower, upper in zip(pvs, lower_pos, upper_pos):
            depth = sig[center, i]
            sig_range = sig[lower:upper+1, i]
            sig_range[np.where(int_sign(sig_range) != int_sign(depth))] = 0.0
            volume = np.sum(sig_range)
            peaks.append(Peak(center=center, lower=lower, upper=upper, depth=depth, volume=volume))
        peaks_list.append(peaks)
    if oned:
        peaks_list = peaks_list[0]
    
    if return_allinfo:
        return peaks_list, cr, cr_crosszero, cr_inflection 
    else:
        return peaks_list

if __name__ == '__main__':
    pass
