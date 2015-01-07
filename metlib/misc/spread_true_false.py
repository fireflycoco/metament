#!/usr/bin/env python

from copy import copy, deepcopy
import numpy as np

__all__ = ['spread_something_1d', 'spread_true', 'spread_false']

def spread_something_1d(arr, something, forward_num=0, backward_num=0):
    """spread some special value in a 1d-array to its neighbourhood.
    """
    res = deepcopy(arr)
    w = np.where(arr == something)
    total_length = arr.shape[0]
    for i in w[0]:
        end = min(total_length, i+forward_num)
        beg = max(0, i-backward_num)
        res[beg:end] = something
    return res

def spread_true(arr, forward_num=0, backward_num=0):
    """spread true value in a array.
    e.g. to filter out neigbours of invalid value:
        bool_arr = value_arr < 0.0
        new_bool_arr = spread_true(bool_arr, forward_num=5, backward_num=1)
        value_arr[np.where(new_bool_arr)] = 0.0
    """
    return spread_something_1d(arr, True, forward_num, backward_num)

def spread_false(arr, forward_num=0, backward_num=0):
    """spread false value in a array.
    e.g. to filter out neigbours of invalid value:
        bool_arr = value_arr > 0.0
        new_bool_arr = spread_false(bool_arr, forward_num=5, backward_num=1)
        valid_value_arr = value_arr[np.where(new_bool_arr)]
    """
    return spread_something_1d(arr, False, forward_num, backward_num)
