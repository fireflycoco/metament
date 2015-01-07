#!/usr/bin/env python
# -*- coding:utf-8 -*-

# lookup_table.py

import os, sys
import numpy as np
import scipy as sp

__all__ = ['lookup_table', 'nc2lut', 'bin2lut', 'enumerate_kwarg_dict', 'strval2lut']

def enumerate_kwarg_dict(**kwargs):
    """Given kwargs of sequences, enumerate_kwarg_dict generates every combination of the kwargs values provided, yielding a dict of the kwargs as keys.

Parameters
----------
kwargs: kwargs of sequences.

Return
------
yields every combination of the kwarg sequences.

Examples
--------
>>> for d in enumerate_kwarg_dict(abc=['a', 'B', 'c'], num=[1, 3]):
...     print d
{'abc': 'a', 'num': 1}
{'abc': 'a', 'num': 3}
{'abc': 'B', 'num': 1}
{'abc': 'B', 'num': 3}
{'abc': 'c', 'num': 1}
{'abc': 'c', 'num': 3}
"""
    keys = []
    lens = []
    for key in kwargs:
        keys.append(key)
        lens.append(len(kwargs[key]))
    for ndind in np.ndindex(*lens):
        nowdict = dict()
        for ind, key in zip(ndind, keys):
            nowdict[key] = kwargs[key][ind]
        yield nowdict

class lookup_table(object):
    def __init__(self, arr, dims):
        """arr is the big ndarray of lookup_table,
        dims is a list of name, values tuple:
        [('dim1',[0.1,0.2,0.3,0.4,0.5]), ('dim2',[35.,40.,50.]),...]
        any dim name is OK.
        """
        self.arr = arr
        # # check
        shape_checker = []
        for d in dims:
            shape_checker.append(len(d[1]))
            l = list(d[1])
            order_checker = sorted(l)
            if l != order_checker and l != reversed(order_checker):
                raise ValueError( \
"""Lookup Table's dim: %s is not in ascending or descending order
    %s
""" % (d[0], d[1]) )
        if tuple(shape_checker) != self.arr.shape:
            raise ValueError("Lookup Table's array and dim mismatch: %s vs %s" \
                    % (self.arr.shape, tuple(shape_checker)))

        self.dims = []
        for d in dims:
            self.dims.append((d[0], np.array(d[1])))
        self.dimnames = [item[0] for item in dims]

    def __call__(self, coordinate=None, **kwargs):
        return self.lookup(coordinate, **kwargs)

    def lookup(self, coordinate=None, **kwargs):
        """lt.lookup(dim1=xxx, dim2=xxx,..., dimn=xxx) or 
        lt.lookup([x, y, z ...])
        For the first method, if len(dims) == len(kwargs), return a scaler value, else, return a lookup_table with len(dims) - len(kwargs) dims.
        For the second method, use the default dim order.
        """
        if coordinate is not None:
            if len(coordinate) > len(self.dimnames):
                raise ValueError("""Too many dims in the given coordinate.
        Lookup Table's dims: %s
        Given coordinate: %s""" % (self.dimnames, coordinate))
            kwargs = dict()
            for i in range(len(coordinate)):
                kwargs[self.dimnames[i]] = coordinate[i]
        
        order = []
        for name in self.dimnames:
            if name in kwargs:
                order.append(name)
        if len(order) != len(kwargs):
            raise ValueError("""Dim names not match.
    Lookup Table's dim names: %s
    Function args: %s""" % ( self.dimnames, kwargs.keys() ) )

        # # get enough good info
        good_info = dict()
        for now_dimname in order:
            # find dim's index in all dim
            for i, dim in enumerate(self.dims):
                if now_dimname == dim[0]:
                    dim_i = i
                    dim_value = np.array(dim[1])
                    break
            now_value = kwargs[now_dimname]
            # # To prevent extrapolating
            if now_value < dim_value.min():
                now_value = dim_value.min()
            elif now_value > dim_value.max():
                now_value = dim_value.max()

            # find layer's index in the dim
            for i in range(len(dim_value)-1):
                if (dim_value[i]<=now_value and now_value<=dim_value[i+1]) \
                    or \
                    (dim_value[i]>=now_value and now_value>=dim_value[i+1]):
                    left_i, right_i = i, i+1
                    break
            left_value, right_value = dim_value[i], dim_value[i+1]
            left_ratio = (right_value - now_value) / (right_value - left_value)
            right_ratio = (now_value - left_value) / (right_value - left_value)
            good_info[now_dimname] = {'dimname':now_dimname, 'dim_i':dim_i, \
                    'left_i':left_i, 'right_i':right_i, \
                    'left_ratio':left_ratio, 'right_ratio':right_ratio }
        # # make a small array for interp
        magic = 'self.arr['
        for dimname in self.dimnames:
            if dimname not in good_info:
                magic += ':,'
            else:
                info = good_info[dimname]
                magic += '%d:%d,' % (info['left_i'], info['right_i']+1)
        magic = magic.rstrip(',') + ']'
#        print magic
        small_arr = eval(magic)
#        print order
        shift = 0
        for dimname in order:
            info = good_info[dimname]
            pos_left_str=[':'] * len(small_arr.shape)
            pos_right_str=[':'] * len(small_arr.shape)
            pos_left_str[info['dim_i']-shift] = '0'
            pos_right_str[info['dim_i']-shift] = '1'
            magic = "small_arr[%s] * info['left_ratio'] + small_arr[%s] * info['right_ratio']" % ( ','.join(pos_left_str), ','.join(pos_right_str) )
#            print magic
            small_arr = eval(magic)
            shift += 1

        if len(order) == len(self.dimnames):
            # in this case, small_arr is already a scaler
            return small_arr
        else:
            # make a smaller lookup_table
            new_dims = []
            for d in self.dims:
                if d[0] not in order:
                    new_dims.append(d)
#            print self.arr.shape
#            print small_arr.shape
            return lookup_table(small_arr, new_dims)

    def fast_lookup(self, pos):
        """lt.lookup(pos)
        pos is a seq of position value in each dim.
        this function does NOT assert whether pos is good.
        And assumes that all dims are in acsending order.
        However, it's not much faster than lookup yet.
        """
        dimlen = len(self.dims)
        left_i = np.zeros(dimlen, dtype='i4')
        right_i = np.zeros(dimlen, dtype='i4')
        left_ratio = np.zeros(dimlen, dtype='f4')
        right_ratio = np.zeros(dimlen, dtype='f4')
        for d_i in range(dimlen):
            now_value = pos[d_i]
            dim_value = self.dims[d_i][1]
            # find layer's index in the dim
            for i in range(len(dim_value)-1):
                if dim_value[i]<=now_value and now_value<=dim_value[i+1]:
                    left_i[d_i], right_i[d_i] = i, i+1
                    break
            left_value, right_value = dim_value[i], dim_value[i+1]
            left_ratio[d_i] = (right_value - now_value) / (right_value - left_value)
            right_ratio[d_i] = 1.0 - left_ratio[d_i]
        
        # # make a small array for interp
        # # Time 100 / 430 mu s
        magic = 'self.arr['
        for d_i in range(dimlen):
            magic += '%d:%d,' % (left_i[d_i], right_i[d_i]+1)
        magic = magic.rstrip(',') + ']'
        small_arr = eval(magic)
# #        small_arr = self.arr
        # # Time 70 / 430 mu s
        for d_i in range(dimlen-1):
            small_arr = small_arr[0,:] * left_ratio[d_i] + \
                        small_arr[1,:] * right_ratio[d_i]
        small_arr = small_arr[0] * left_ratio[d_i] + \
                        small_arr[1] * right_ratio[d_i]
        # # Time
        return small_arr

    def reverse_lookup(self, value, span_num=11):
        """find corresponding dim position for the value.
    only for 1d lookup table.
    when 2 adjacent values are identical, 
    return span_num values between the 2 dim values"""
        if len(self.dims) != 1:
            raise RuntimeError("lookup_table.reverse_lookup(value) only \
works for 1-d lookup_tables. This one is %d-d." % len(self.dims))
        if span_num <= 1:
            span_num = 2
        res = []
        diffs = self.arr[1:] - self.arr[:-1]
        for i in range(len(self.arr)-1):
            if value >= np.min(self.arr[i:i+2]) and \
                    value <= np.max(self.arr[i:i+2]):
                try:
                    left_ratio = (self.arr[i+1] - value) / diffs[i]
                    right_ratio = (value - self.arr[i]) / diffs[i]
                    dim_value = self.dims[0][1][i] * left_ratio + \
                            self.dims[0][1][i+1] * right_ratio
                    if np.isfinite(dim_value):
                        res.append(dim_value)
                    else: 
                        raise ZeroDivisionError
                except ZeroDivisionError:
                    res.extend(np.linspace(self.dims[0][1][i],
                        self.dims[0][1][i+1],span_num))
        return np.unique(np.array(res))

    def __str__(self):
        res = 'lookup_table'
        for dname, diminfo in self.dims:
            res = res + '\n%8s : %s' % (dname, diminfo)
        return res

def nc2lut(ncfname):
    pass

def bin2lut(descfname, **kwargs):
    """Load a raw binary file as lut directly. 
    Parameters:
        descfname: a description txt file in the format below:
            # Comment lines starts with #
            binary_file : /path/to/binary_file
            dtype : datatype of the lut, a character in 'cb1silfdFD'
            dimname1 : 1 2 3 4 5
            dimname2 : 10.0, 20.0, 30.0, 40.0
            ...
            # dims can either be seperated with ',' or ' '
        kwargs: 
            mem_type:
                a character in 'cb1silfdFD'
            byteswap:
                0 : no swap
                1 : swap
"""
    from scipy.io.numpyio import fread
    descf_folder_path = os.path.dirname(descfname)
    descf_basename = os.path.basename(descfname)
    descf = open(descfname)
    dims = []
    shp = []
    binary_file = '%s.bin' % descf_basename
    dtype = 'f'
    for l in descf:
        l = l.strip()
        if len(l) == 0 or l.startswith('#'):
            continue
        name, value = l.split(':', 1)
        name = name.strip()
        value = value.strip()
        if name == 'binary_file':
            binary_file = value
        elif name == 'dtype':
            dtype = value
        else:
            if ',' in value:
                sep = ','
            else:
                sep = ' '
            arrayvalue = np.fromstring(value, sep=sep, dtype='f4')
            dims.append((name, arrayvalue))
            shp.append(arrayvalue.shape[0])
    total_elements_num = np.multiply.reduce(shp)
    shp = tuple(shp)
    binary_file_abspath = os.path.join(descf_folder_path, binary_file)
    f = open(binary_file_abspath, 'rb')
    lut_arr = fread(f, total_elements_num, dtype).reshape(shp)
    f.close()
    lut = lookup_table(lut_arr, dims)
    return lut

def strval2lut(dimstrs, values, sep='_', func=None):
    """strval2lut builds LUT with a seq of str which contains dimension information and a seq of values.

    Parameters
    ----------
    dimstrs: a seq of str in the format of dim1=val1_dim2=val2_dim3=val3... , where '_' is used as default separater. This seq should be sorted.
    values:  corresponding value seq .
    sep: separater, default is '_' (in this case, do not use '_' in dim names).
    func: an optional func to process the values, default is None. 

    Return
    ------
    A lookup_table.

    Examples
    --------
    >>> dimstrs = ['abc=-1.0_num=1', 'abc=-1.0_num=3', 'abc=0.0_num=1', 'abc=0.0_num=3', 'abc=1.5_num=1', 'abc=1.5_num=3']
    ... lut = strval2lut(dimstrs, [0,1,2,3,4,5])
    ... print lut
    ... print lut.arr
    lookup_table
         abc : [-1.   0.   1.5]
         num : [ 1.  3.]
    [[0 1]
     [2 3]
     [4 5]]
    """
    if isinstance(func, np.ufunc):
        true_values = func(values)
    elif callable(func):
        true_values = [func(item) for item in values]
    else:
        true_values = values
    true_values = np.array(true_values)
    dim_eq_val_tokens = dimstrs[0].split(sep)
    dimnames = [d_e_v.split('=')[0] for d_e_v in dim_eq_val_tokens]
    dimdict = dict()
    dimsizes = []
    diminfos = []
    for dimname in dimnames:
        dimdict[dimname] = set()
    for dimstr in dimstrs:
        dvs = dimstr.split(sep)
        for dv in dvs:
            dim, val = dv.split('=')
            val = float(val)
            dimdict[dim].add(val)
    for dimname in dimnames:
        diminfos.append((dimname, sorted(dimdict[dimname])))
        dimsizes.append(len(dimdict[dimname]))
    finalshape = tuple(dimsizes)
    bigarr = true_values.reshape(finalshape)
    lut = lookup_table(bigarr, diminfos)
    return lut

if __name__ == '__main__':
    # # Test code
#    dim1 = np.array([0.1,0.2,0.3,0.4,0.5])
#    dim2 = np.array([35,40,50])
#    lt_arr = np.arange(len(dim1)*len(dim2)).reshape((len(dim1),len(dim2)))
#    lt = lookup_table(lt_arr, [('dim1',dim1), ('dim2',dim2)])
#    print dim1, dim2
#    print lt_arr
#    res1 = lt.lookup(dim1=0.45, dim2=45.0)
#    res1_fast = lt.fast_lookup((0.45, 45.0))
#    print res1, res1_fast
#
#    res2 = lt.lookup(dim1=0.35)  # res2 is a smaller lookup_table
#    print res2.arr
#    print res2.dims
#    print res2.lookup(dim2=38.0)
#    print res2.reverse_lookup(8.1)
#    print res2.reverse_lookup(8.2)
#    print res2.reverse_lookup(8.3)
#    # # The follow line will trigger an exception 
#    # # because there's no dimension named 'dim3'
#    res3 = res2.lookup(dim3=45.0)
#    print res3
    
    dstrs = []
    for d in enumerate_kwarg_dict(abc=[-1.0, 0.0, 1.5], num=[1, 3]):
        print d
        dstrs.append('%s=%s_%s=%s' % ('abc', d['abc'], 'num', d['num']))
    print dstrs
    lut = strval2lut(dstrs, range(len(dstrs)))
    print lut
    print lut.arr
