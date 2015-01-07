#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np
import re
from matplotlib import mlab

__all__ = ['load_aeronet']

def load_aeronet(fname, keep_fields='all', header=False):
    """loads aeronet lev 2.0 csv file.
    fname: data file name
    keep_fields: 'all' or a list of fields
    header: whether to return header information along with the data.
    """
    std_day = datetime(1900,1,1,0,0,0)
    def date2daynum(datestr):
        the_day = datetime.strptime(datestr, '%d:%m:%Y')
        return float((the_day - std_day).days)

    def time2seconds(timestr):
        h, m, s = [int(t) for t in timestr.split(':')]
        return float(h * 3600 + m * 60 + s)

    def daynum_seconds2datetime(daynum, seconds):
        return std_day + timedelta(days=int(daynum), seconds=int(seconds))

    headlines = []
    f = open(fname, 'r')
    for line_i, line in enumerate(f):
        line = line.rstrip()
        if line.startswith('Date(dd-mm-yy'):
            datefield, timefield = [re.sub(r'\W', '', tk) for tk in line.split(',')[0:2]]
            break
        headlines.append(line)
    skip_header_lines = line_i

    if header:
        headline = ','.join(headlines)
        headerd = dict()
        for attrname, converter in [('location', str), ('long', float), ('lat', float), ('elev', float), ('nmeas', int), ('PI', str), ('email', str)]:
            m = re.search(r'%s.{0,1}=([^,\s]*)' % attrname, headline, flags=re.I)
            if m:
                try:
                    headerd[attrname] = converter(m.group(1))
                except Exception:
                    pass

    rawd = np.genfromtxt(fname, skip_header=skip_header_lines, delimiter=',', names=True, converters={0:date2daynum, 1:time2seconds})
    lend = len(rawd)
    dates = np.zeros(len(rawd), dtype='O')
    for i in range(lend):
        dates[i] = daynum_seconds2datetime(rawd[datefield][i], rawd[timefield][i])

    newd = mlab.rec_append_fields(rawd, 'datetime', dates)
    newd = mlab.rec_drop_fields(newd, [datefield, timefield, 'Last_Processing_Date'])

    if keep_fields is not 'all':
        keep_fields = ['datetime'] + keep_fields
#        print keep_fields
        newd = mlab.rec_keep_fields(newd, keep_fields)
    if header:
        return newd, headerd
    else:
        return newd
