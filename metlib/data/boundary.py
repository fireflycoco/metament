#!/usr/bin/env python2.7

# boundary.py

import os, sys
import numpy as np
from datetime import timedelta
from metlib.misc.datatype import isinteger

__all__ = ['Limiter', 'intersection']

class Limiter(object):
    def __init__(self, min_value, max_value):
        self.min_value = min_value
        self.max_value = max_value

    def __call__(self, orig):
        if isinstance(orig, np.ndarray):
            result = orig.copy()
            if self.min_value is not None:
                result[result < self.min_value] = self.min_value
            if self.max_value is not None:
                result[result > self.max_value] = self.max_value
        elif isinstance(orig, list):
            result = np.array(orig)
            if self.min_value is not None:
                result[result < self.min_value] = self.min_value
            if self.max_value is not None:
                result[result > self.max_value] = self.max_value
            result = list(result)
        elif isinstance(orig, slice):
            step  = 1 if orig.step is None else orig.step
            start = 0 if orig.start is None else orig.start
            if self.min_value is not None and start < self.min_value:
                start = self.min_value + (step - (self.min_value - start) % step) % step
            if self.max_value is not None and start > self.max_value:
                return slice(self.max_value, self.max_value, 1)
            stop = orig.stop
            if stop is not None and self.max_value is not None and stop > self.max_value+1:
                stop = self.max_value+1
            result = slice(start, stop, step)
        else:
            if self.min_value is not None and orig < self.min_value:
                result = self.min_value
            elif self.max_value is not None and orig > self.max_value:
                result = self.max_value
            else:
                result = orig
        return result

    def __getitem__(self, orig_slice):
        return self.__call__(orig_slice)

    def __repr__(self):
        return "Limiter[%s,%s]" % (self.min_value, self.max_value)
    
def intersection(a, b):
    """calculate the length of 2 lines' intersection.
    a, b: (beg, end) tuple of two lines; or two arrays.
    returns: intersection length of the two lines.
    e.g., intersection( (3, 6), (4, 8) ) -> 2, which is 4 ~ 6.
    """
    ab = np.r_[a, b]
    res = np.abs(a[-1] - a[0]) + np.abs(b[-1] - b[0]) - (np.max(ab) - np.min(ab))
    if isinstance(res, timedelta):
        zerotd = timedelta(seconds=0)
        if res < zerotd:
            res = zerotd
    elif isinteger(res):
        if res < 0:
            res = 0
    else:
        if res < 0.0:
            res = 0.0
    return res

if __name__ == '__main__':
    pass
