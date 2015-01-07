#!/usr/bin/env python

import numpy as np
import sys
from .logical import *
from .misc import *
from .parser import *

__all__ = ['datetime_split', 
        'split_type',
        'split_season', 'split_month', 
        'split_weekday', 'split_hour', 'split_year',
        'split_year_season', 'split_year_month',
        'split_datetime_str'
        ]

def datetime_split(rec, dts=None, funcs=[]):
    if dts is None:
        dts = rec['datetime']
    dts = parse_datetime(dts)
    bins = []
    for func in funcs:
        part = rec[np.where(func(dts))]
        bins.append(part)
    return bins

def split_type(rec, codes, return_info=False):
    unique_codes = list(set(codes))
    funcs = [lambda s, code=code: s == code for code in unique_codes]
    res = datetime_split(rec, np.array(codes), funcs)
    if return_info:
        return res, unique_codes
    else:
        return res

def split_season(rec, dts=None):
    return datetime_split(rec, dts, funcs=[
        lambda s: season_is('spring', s),
        lambda s: season_is('summer', s),
        lambda s: season_is('autumn', s),
        lambda s: season_is('winter', s)
        ]
        )

def split_month(rec, dts=None):
    return datetime_split(rec, dts, 
            funcs=[lambda s, month=month: month_is(month, s) for month in range(1, 13)])

def split_weekday(rec, dts=None):
    return datetime_split(rec, dts, 
            funcs=[lambda s, weekday=weekday: weekday_is(weekday, s) for weekday in range(1, 8)] )

def split_hour(rec, dts=None):
    return datetime_split(rec, dts, 
            funcs=[lambda s, hour=hour: hour_is(hour, s) for hour in range(0,24)] )

def split_year(rec, dts=None, start_year=None, end_year=None, return_info=False):
    if dts is None:
        dts = rec['datetime']
    if start_year is None:
        start_year = dts[0].year
    if end_year is None:
        end_year = dts[-1].year + 1
    years = range(start_year, end_year)
    res = datetime_split(rec, dts, 
            funcs=[lambda s, year=year: year_is(year, s) for year in years] )
    if return_info:
        return res, years
    else:
        return res

def split_year_season(rec, dts=None, start_year=None, end_year=None, start_season=None, end_season=None, return_info=False):
    if dts is None:
        dts = rec['datetime']
    auto_start_year = False
    auto_end_year = False
    if start_year is None:
        start_year = dts[0].year
        auto_start_year = True
    if end_year is None:
        end_year = dts[-1].year + 1
        auto_end_year = True
    if start_season is None:
        start_season = month2season(dts[0].month, '1234')
    if end_season is None:
        end_season = month2season(dts[-1].month, '1234')
    if start_season == 4  and auto_start_year:
        start_year -= 1
    if end_season == 1 and auto_end_year:
        end_year -= 1
    ys = [(y, s) for y in range(start_year, end_year) for s in range(1,5)]
    end_i = -(4-end_season)
    if end_i == 0: end_i = None 
    ys = ys[start_season-1:end_i]
    res = datetime_split(rec, dts, 
            funcs=[lambda s, year=year, season=season: year_season_is(year, season, s) for year, season in ys] )
    if return_info:
        return res, ys
    else:
        return res

def split_year_month(rec, dts=None, start_year=None, end_year=None, start_month=None, end_month=None, return_info=False):
    if dts is None:
        dts = rec['datetime']
    if start_year is None:
        start_year = dts[0].year
    if end_year is None:
        end_year = dts[-1].year + 1
    if start_month is None:
        start_month = dts[0].month
    if end_month is None:
        end_month = dts[-1].month
    ym = [(y, m) for y in range(start_year, end_year) for m in range(1,13)]
    end_i = -(12-end_month)
    if end_i == 0: end_i = None
    ym = ym[start_month-1:end_i]
    res = datetime_split(rec, dts, 
            funcs=[lambda s, year=year, month=month: year_is(year, s) & month_is(month, s) for year, month in ym] )
    if return_info:
        return res, ym
    else:
        return res

def split_datetime_str(rec, dts=None, fmt='%Y%m%d', return_method='index'):
    """split [rec]array by datetime_fmt.
    dts: if not None, use this array as datetime instead of rec's own 'datetime' field.
    return_method: 'index' | 'item' . 
    """
    if dts is None:
        dts = rec['datetime']

    dtstrd = {}
    for i, dt in enumerate(dts):
        dtstr = dt.strftime(fmt)
        pool = dtstrd.get(dtstr, None)
        if pool is None:
            dtstrd[dtstr] = [i]
        else:
            pool.append(i)
    
    if return_method == 'index':
        return dtstrd
    elif return_method == 'item':
        rec = np.array(rec)
        res = {}
        for dtstr, indice in dtstrd.iteritems():
            res[dtstr] = rec[indice]
        return res

