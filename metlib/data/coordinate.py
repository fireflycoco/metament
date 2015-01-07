#!/usr/bin/env python

# coordinate.py
"""This module provides functions on transforming coordinate"""

import os, sys
import numpy as np

__all__ = ['CoordTransformer1D', 'PointTransformer'] 

class CoordTransformer1D(object):
    def __init__(self, old_x1, old_x2, new_x1, new_x2):
        self.old_x1 = old_x1
        self.old_x2 = old_x2
        self.new_x1 = new_x1
        self.new_x2 = new_x2
        self._calc_scale()

    def _calc_scale(self):
        self.scale = (self.new_x2 - self.new_x1) / (self.old_x2 - self.old_x1)
        self.offset =  self.new_x1 - self.scale * (self.old_x1) 

    def __call__(self, x, reverse=False):
        if reverse:
            return (x - self.offset) / self.scale
        else:
            return x * self.scale + self.offset

class PointTransformer(object):
    def __init__(self, old_x1, old_x2, new_x1, new_x2, old_y1, old_y2, new_y1, new_y2):
        self.xtransformer = CoordTransformer1D(old_x1, old_x2, new_x1, new_x2)
        self.ytransformer = CoordTransformer1D(old_y1, old_y2, new_y1, new_y2)

    def __call__(self, x, y, reverse=False):
        return self.xtransformer(x, reverse), self.ytransformer(y, reverse)

if __name__ == '__main__':
    tester = CoordTransformer1D(1.0, 2.0, 1.0, 3.0)
    print tester(0)
    print tester(1)
    print tester(2, reverse=True)
