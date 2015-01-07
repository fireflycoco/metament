#!/usr/bin/env python

from datetime import datetime, timedelta
from dateutil.parser import parse
from .parser import parse_datetime, parse_timedelta

__all__ = ['datetime_range']

def datetime_range(beg, end, tdelta):
    """Returns a list of datetimes from beg to end with tdelta.
    """
    if not isinstance(beg, datetime):
        beg = parse_datetime(beg)
    if not isinstance(end, datetime):
        end = parse_datetime(end)
    if not isinstance(tdelta, timedelta):
        tdelta = parse_timedelta(tdelta)
    result = []
    now = beg
    while now < end:
        result.append(now)
        now = now + tdelta
    return result

if __name__ == "__main__":
    print datetime_range(20110101, 20110105, '1d')
