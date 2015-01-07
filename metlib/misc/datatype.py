#!/usr/bin/env python2.7

# datatype.py

import os, sys
import collections
from array import array as pythonarray
import numpy as np

__all__ = ['limited_int', 'Singleton', 'Null', 'NullClass', 
            'isnull', 'isinteger', 'isfloat', 'isseq', 'isstr']

class limited_int(object):
    """An integer-like class with max and min limits."""
    def __init__(self, value, vmin, vmax):
        self.value = long(value)
        self.set_limit(vmin, vmax)

    def set_value(self, value):
        self.value = long(value)
        self._check()

    def rel_range(self, beg, end, step=1, under="min", over="max"):
        if under == 'min':
            under = self.vmin
        if over == 'max':
            over = self.vmax
        res = np.arange(beg, end, step, dtype='i8') + self.value
        res[res > self.vmax] = over
        res[res < self.vmin] = under
        return res
    
    def set_limit(self, vmin, vmax):
        self.vmin = long(vmin)
        self.vmax = long(vmax)
        self._check()

    def get_limit(self):
        return self.vmin, self.vmax

    def _check(self):
        if self.value > self.vmax:
            self.value = self.vmax
        if self.value < self.vmin:
            self.value = self.vmin

    def __iadd__(self, other):
        self.value += other
        self._check()
        return self

    def __isub__(self, other):
        self.value -= other
        self._check()
        return self

    def __int__(self):
        return int(self.value)

    def __long__(self):
        return long(self.value)

    def __float__(self):
        return float(self.value)

    def __lt__(self, other):
        return int(self) < other

    def __le__(self, other):
        return int(self) <= other

    def __gt__(self, other):
        return int(self) > other

    def __ge__(self, other):
        return int(self) >= other

    def __eq__(self, other):
        return int(self) == other

    def __ne__(self, other):
        return int(self) != other

    def __add__(self, other):
        return int(self) + other

    def __sub__(self, other):
        return int(self) - other

    def __mul__(self, other):
        return int(self) * other

    def __div__(self, other):
        return int(self) / other

    def __radd__(self, other):
        return other + int(self) 

    def __rsub__(self, other):
        return other - int(self)

    def __rmul__(self, other):
        return other * int(self) 

    def __rdiv__(self, other):
        return other / int(self)

    def __str__(self):
        return '%d:<%d %d>' % (self.value, self.vmin, self.vmax)

    def __repr__(self):
        return self.__str__()

class Singleton(object):
    """A do-nothing class. 
    From A. Martelli et al. Python Cookbook. (O'Reilly)
    Thanks to Juergen Hermann."""
    def __new__(cls, *args, **kwargs):
        if '_inst' not in vars(cls):
            cls._inst = super(Singleton, cls).__new__(cls, *args, **kwargs)
        return cls._inst

class NullClass(Singleton):
    """A do-nothing class. 
    From A. Martelli et al. Python Cookbook. (O'Reilly)
    Thanks to Dinu C. Gherman, Holger Krekel.
    """
    def __init__(self, *args, **kwargs): pass
    def __call__(self, *args, **kwargs): return self
    def __repr__(self): return "Null"
    def __nonzero__(self): return False
    def __getattr__(self, name): return self
    def __setattr__(self, name, value): return  self
    def __delattr__(self, name): return self
    def __len__(self): return 0
    def __iter__(self): return iter(())
    def __getitem__(self, i): return self
    def __delitem__(self, i): return self
    def __setitem__(self, i): return self

Null = NullClass()

def isnull(value):
    return isinstance(value, (type(None), NullClass))

def isinteger(value):
    return isinstance(value, (int, long, np.integer))

def isfloat(value):
    return isinstance(value, (float, np.float))

def isseq(value):
    return isinstance(value, (list, tuple, bytearray, buffer, xrange, pythonarray,
        collections.deque, 
        np.ndarray, np.recarray,))

def isstr(value):
    return isinstance(value, (str, unicode))

if __name__ == '__main__':
    a = limited_int(3, 0, 4)
    print a
    print a.rel_range(-4, 5)
    a += 1
    print a
    a += 1 
    print a
    a += 1 
    print a
    print int(a)+3
