#!/usr/bin/env python
# -*- coding:utf-8 -*-

# misc.py

#import os, sys
#import re
from datetime import datetime, timedelta
import numpy as np
#import scipy as sp
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
from matplotlib.dates import date2num as mpl_date2num
from matplotlib.dates import num2date as mpl_num2date
from metlib.misc.datatype import isseq

from .parser import *
__all__ = ['month2season', 'datetime2season', 'datetime2yearseason',
        'str2datetime', 'datetime_match', 'datetime_filter',
        'season_names', 'month_names', 'month_short_names', 
        'matlab_num2date', 'matlab_date2num',
        'mpl_date2num', 'mpl_num2date',
        ]

season_names = ['Spring', 'Summer', 'Autumn', 'Winter']
month_names  = ['January', 'February', 'March', 'April', 'May', 'June',
        'July', 'August', 'September', 'October', 'November', 'December']
month_short_names = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

def month2season(month, outformat='0123'):
    """Converts month number ( 1-12 ) to season number or names.
    Parameters:
        month : int or seq of int ( 1 - 12 )
        outformat: '0123', '1234', 'name' or ANY seq with at least 4 elements
"""
    if outformat == 'name':
        season_names = np.array(['Spring', 'Summer', 'Autumn', 'Winter'], dtype='O')
    elif outformat == '0123':
        season_names = np.array([0,1,2,3], dtype='i4')
    elif outformat == '1234':
        season_names = np.array([1,2,3,4], dtype='i4')
    else:
        if len(outformat) < 4:
            raise ValueError, 'outformat: "%s" length less than 4' % outformat
        season_names = np.array(list(outformat))

    if isseq(month):
        season_index = (np.array(month, dtype=int)+9) / 3 % 4
    else:
        season_index = (int(month) + 9) / 3 % 4
    return season_names[season_index]

def datetime2season(dts, outformat='0123'):
    """Converts datetime to season number or names.
    Parameters:
        dts: datetime / datetime seq.
        outformat: '0123', '1234', 'name' or ANY seq with at least 4 elements
"""
    if isseq(dts):
        months = [dt.month for dt in T(dts)]
    else:
        months = T(dts).month
    return month2season(months, outformat=outformat)

def datetime2yearseason(dts, seasonformat='name', sep='_', DJF='JF_year'):
    """Converts datetime to "Year_Season" strs
    Parameters:
        dts: datetime / datetime seq;
        seasonformat: '0123', '1234', 'name' or ANY seq with at least 4 elements;
        sep: result str is "Year" + sep + "Season";
        DJF: 
            'JF_year' (default): December belongs to the next year (same as Jan/Feb).
            'D_year': Jan/Feb belongs to the prev year (same as Dec).
            '': Dec/Jan/Feb belongs to its own year.
    Returns:
        "Year_Season" str/seq.
"""
    dts = T(dts)
    if isseq(dts):
        scalar_res = False
    else:
        scalar_res = True
        dts = np.array([dts])

    seasons = datetime2season(dts, seasonformat)
    years = np.array([dt.year for dt in dts])
    months = np.array([dt.month for dt in dts])
    if DJF == 'JF_year':
        years[months == 12] += 1
    elif DJF == 'D_year':
        years[months == 1] -= 1
        years[months == 2] -= 1
    res = np.array(['%s%s%s' % (year, sep, season) for year, season in zip(years, seasons)])
    if scalar_res:
        res = res[0]
    return res

def str2datetime(fmt='%Y-%m-%d %H:%M:%S', *arrays ):
    """Converts str arrays into datetime arrays.
    Parameters:
        fmt: datetime format string. When using more than one array as input, extra white space should be added. 
        arrays: string or int arrays to parse.
    Returns:
        datetime array
    """
    # TODO: int array not working well now. add something to format int arrays into str first.
    n = len(arrays)
    l = len(arrays[0])
    res = np.zeros(l, dtype='O')
    for i in range(l):
        try:
            datestrs = []
            for j in range(n):
                datestrs.append(str(arrays[j][i]))
            datestr = ' '.join(datestrs)
            res[i] = datetime.strptime(datestr, fmt)
        except Exception, e:
            print e
            res[i] = None
    return res

def datetime_match(rec, ref_dts, fmt="%Y%m%d%H%M%S", fmt2=None, rec_dts_field='datetime', return_index=False, fill_value=np.nan):
    """match a subset of rec to match ref_dts
    rec: recarray to select from.
    ref_dts: reference datetime array.
    fmt: common datetime string format for matching.
    fmt2: ref_dts's format string. if None, use fmt
    rec_dts_field: 'datetime', 'date', etc. Use None if rec is a datetime seq. Use a datetime seq if rec does not contain a dts field.
    return_index: if False: return selected rec only; if True: also return matching index.
    fill_value: fill value.
Bug:  There may be None in returned index.
    """
    if fmt2 is None:
        fmt2 = fmt
    rec = np.array(rec)
    ref_dts = parse_datetime(ref_dts)
    if rec_dts_field is None:
        rec_dts = rec
    elif isinstance(rec_dts_field, (str, unicode)):
        rec_dts = rec[rec_dts_field]
    else:
        rec_dts = rec_dts_field
    rec_dtstrs = [dt.strftime(fmt) for dt in rec_dts]
    rec_dict = dict(zip(rec_dtstrs, range(len(rec_dtstrs))))
    ref_dtstr = [dt.strftime(fmt2) for dt in ref_dts]
    res_i = []
    for dtstr in ref_dtstr:
        i = rec_dict.get(dtstr, None)
        res_i.append(i)
    if rec_dts_field is None:
        res = np.zeros(len(res_i), dtype='O')
    else:
        res = np.atleast_1d(np.zeros(len(res_i), rec.dtype))
        fields = rec.dtype.names
    res[:] = fill_value
    if rec_dts_field is None:
        for i, rec_i in enumerate(res_i):
            if rec_i is None:
                pass
            else:
                res[i] = rec[rec_i]
    else:
        for i, rec_i in enumerate(res_i):
            if fields is None:
                if rec_i is None:
                    res[i] = fill_value
                else:
                    res[i] = rec[rec_i]
            else:
                if rec_i is None:
                    # TODO right hand stuff still not perfect
                    if isinstance(rec_dts_field, (str, unicode)):
                        res[rec_dts_field][i] = parse_datetime(ref_dtstr[i])
                else:
                    for f in fields:
                        res[f][i] = rec[f][rec_i]
    if return_index:
        return res, np.array(res_i)
    else:
        return res

def datetime_filter(rec, ref_dts, fmt="%Y%m%d%H%M%S", fmt2=None, rec_dts_field='datetime', return_index=False, fill_value=np.nan):
    """filter a subset of rec to match ref_dts
    rec: recarray to select from.
    ref_dts: reference datetime array. if it's a seq of str AND fmt2 is None: use it directly
    fmt: common datetime string format for matching.
    fmt2: ref_dts's format string. if None, use fmt
    rec_dts_field: 'datetime', 'date', etc. Use None if rec is a datetime seq. Use a datetime seq if rec does not contain a dts field.
    return_index: if False: return selected rec only; if True: also return matching index.
    fill_value: fill value.
    """
    rec = np.array(rec)
    if isinstance(ref_dts[0], (str, unicode)) and fmt2 is None:
        ref_dtstr = ref_dts
    else:
        if fmt2 is None:
            fmt2 = fmt
        ref_dts = parse_datetime(ref_dts)
        ref_dtstr = [dt.strftime(fmt2) for dt in ref_dts]
    if rec_dts_field is None:
        rec_dts = rec
    elif isinstance(rec_dts_field, (str, unicode)):
        rec_dts = rec[rec_dts_field]
    else:
        rec_dts = rec_dts_field
    rec_dtstrs = [dt.strftime(fmt) for dt in rec_dts]
    rec_dict = dict()
    for i in range(len(rec_dtstrs)):
        if rec_dtstrs[i] not in rec_dict:
            rec_dict[rec_dtstrs[i]] = [i]
        else:
            rec_dict[rec_dtstrs[i]].append(i)
    res_i = []
    for dtstr in ref_dtstr:
        i_list = rec_dict.get(dtstr, None)
        if i_list is not None:
            res_i.extend(i_list)
    res = np.atleast_1d(rec[(res_i,)])
    if return_index:
        return res, np.array(res_i)
    else:
        return res

def matlab_date2num(the_dates):
    """
    """
    return mpl_date2num(parse_datetime(the_dates)) + 366.0

def matlab_num2date(the_nums, **kwargs):
    """
    """
    if np.isscalar(the_nums):
        return mpl_num2date(the_nums - 366.0, **kwargs)
    else:
        return mpl_num2date(np.array(the_nums) - 366.0, **kwargs) 
