#!/usr/bin/env python
import sys
from datetime import datetime, date

__all__ = ['rec2csv_neat']

def rec2csv_neat(rec, f, formatd={}, delimiter=','):
    """rec2csv_neat is a neat version of mlab.rec2csv.
    rec: recarray to save.
    formatd: dict of formats for each fields, who's values may be format str, function, etc. For field 'datetime','date','time', '%Y%m%d...' style format can be used"
    """
    dtp = rec.dtype
    fieldnames = dtp.names
    formaters = {}
    for fn in fieldnames:
        if fn in formatd:
            the_fmt = formatd[fn]
            if isinstance(the_fmt, (str, unicode)):
                if isinstance(rec[fn][0], (date, datetime)) :
                    formaters[fn] = lambda dt: dt.strftime(formatd[fn])
                else:
                    formaters[fn] = lambda val: formatd[fn] % (val, )
            elif callable(the_fmt):
                formaters[fn] = the_fmt
            else:
                sys.stderr.write('Warning: given formater for field %s : %s  is not callable. Not using.\n' % (fn, the_fmt))
                formaters[fn] = str
        else:
            formaters[fn] = str
    if isinstance(f, (str, unicode)):
        f = open(f, 'w')
    f.write(delimiter.join(fieldnames))
    f.write('\n')
    for i in range(len(rec)):
        try:
            to_write = [formaters[fn](rec[fn][i]) for fn in fieldnames]
            final_str = delimiter.join(to_write)
            f.write(final_str)
            f.write('\n')
                
        except Exception as e:
            sys.stderr.write('%s\n' % e)

#if __name__ == "__main__":
#    from matplotlib import mlab
#    rec = mlab.csv2rec('./test.csv')
#    print rec
#    rec2csv_neat(rec, './testout.csv', formatd={'datetime':"%Y%m%d", 'lc':'%.3f','std':'%.9f'})
