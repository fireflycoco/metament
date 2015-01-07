import numpy as np
from netCDF4 import Dataset

__all__ = ['load_nc_profile']

def load_nc_profile(fname, length=None, fill_value=0.0):
    f = Dataset(fname)
    rawdata = f.variables['data'][:]
    orig_length  = rawdata.shape[-1]
    if length is not None:
        new_shape = list(rawdata.shape[:-1])
        new_shape.append(length)
        result = np.zeros(new_shape)
        if length > orig_length:
            result[...,:orig_length] = rawdata
            result[...,orig_length:] = fill_value
        else:
            result[...,:length] = rawdata[...,:length]
        return result
    else:
        return rawdata
