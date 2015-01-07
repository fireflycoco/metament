#!/usr/bin/env python

import re
from datetime import datetime, timedelta, date
from dateutil.relativedelta import relativedelta
from dateutil.parser import parse
import numpy as np

__all__ = ['parse_datetime', 'parse_timedelta', 'T', 'TD']

def parse_datetime(timestr, force_datetime=True):
    """Try parse timestr or integer or list of str/integer into datetimes
    timestr: single value or seq of:
                int/str: YYYYMMDD[HH[MM[SS]]] , delimiters are also allowed, 
                    e.g.: YYYY-MM-DD HH:MM:SS.
                    Julian Days, e.g., YYYYJJJ[...] is also OK.
                datetime/date
    force_datetime: 
                if True: return datetime even if input is date.
                if False: return date if input is date.
    """
    numonly_pattern = re.compile(r'\D*')
    jfmtstr_all = '%Y%j%H%M%S'
    gfmtstr_all = '%Y%m%d%H%M%S'
    return_single = False
    if isinstance(timestr, (int, long, np.integer, str, unicode, datetime, date)):
        timestrs = [timestr]
        return_single = True
    else:
        timestrs = timestr
    res_list = []
    for timestr in timestrs:
        if isinstance(timestr, (datetime, )):
            res_list.append(timestr)
            continue
        if isinstance(timestr, (date, )):
            if force_datetime == True:
                res_list.append(datetime.fromordinal(timestr.toordinal()))
            else:
                res_list.append(timestr)
            continue
        if isinstance(timestr, (int, long, np.integer)):
            timestr = str(timestr)
            
        if len(timestr) == 10:
            if re.match(r'^\d{10}$', timestr):
                timestr = timestr + '00'
            elif re.match(r'^\d{4}.\d{2}.\d{2}$', timestr):
                timestr = timestr + ' 00:00:00'
        try:
            res = parse(timestr)
        except ValueError:
            try:
                numonly_str = re.sub(numonly_pattern, '', timestr)
                nl = len(numonly_str)
                if nl >= 7:
                    if nl % 2 == 1:   #YYYYJJJ[...]
                        fmtstr = jfmtstr_all[:4+nl-7]
                    else:             #YYYYMMDD[...]
                        fmtstr = gfmtstr_all[:4+nl-8]
                    res = datetime.strptime(numonly_str, fmtstr)
                else:
                    res = parse(numonly_str)
            except Exception:
                res = None
        res_list.append(res)
    if return_single:
        return res_list[0]
    else:
        return res_list

_tdelta_dict = {'d':'days', 'h':'hours', 'm':'minutes', 's':'seconds'}
_rel_tdelta_dict = {'Y':'years', 'M':'months', 'W':'weeks'}
def parse_timedelta(timestr):
    """Try parse single value or seq of timestr/integer/timedelta into timedeltas.
    e.g. parse_timedelta(['30m', '2h', '1d', '15s']).
    Y: years
    M: months
    W: weeks
    d: days
    h: hours
    m: minutes
    s: seconds (default)
    """
    return_single = False
    if isinstance(timestr, (int, long, np.integer, str, unicode, timedelta)):
        timestrs = [timestr]
        return_single = True
    else:
        timestrs = timestr
    res_list = []
    for timestr in timestrs:
        if isinstance(timestr, timedelta):
            res_list.append(timestr)
            continue
        if isinstance(timestr, (int, long, np.integer)):
            timestr = str(timestr)
        timestr = timestr.strip()
        m = re.match(r'([-+0-9]+)([dhmsYMW]*)', timestr)
        if not m:
            res = None
        else:
            value = int(m.group(1))
            symbol = m.group(2)
            if symbol in _rel_tdelta_dict:
                res = relativedelta(**{_rel_tdelta_dict[symbol]:value})
            else:
                if symbol == '':
                    symbol = 's'
                res = timedelta(**{_tdelta_dict[symbol]:value})
        res_list.append(res)
    if return_single:
        return res_list[0]
    else:
        return res_list

T = parse_datetime
TD = parse_timedelta

if __name__ == '__main__':
    print parse_datetime(20090101)
    print parse_datetime([20110101, 2011010201, 20110103010101])
    print parse_datetime(['20110101', '2011010201', '20110103010101'])
    print parse_timedelta('35d')
    print parse_timedelta([35, 40, 2, -3])
    print parse_timedelta(['35', '40m', '2d', '-3h'])
