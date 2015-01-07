from netCDF4 import Dataset
from pyproj import Proj
import numpy as np

"""Warning: this module is still buggy. 
"""

class M3_Data_Pos(object):
    def __init__(self):
        self.i = 0
        self.j = 0
        self.k = 0
        self._points = None
        self._weights = None


class M3_File(Dataset):
    """M3_File contains several useful attrs and functions.
    attrs:
        f: a netCDF4.Dataset object which contains everything in the M3 file
        proj: a pyproj.Proj object for coordinate conversion
    """

    def __init__(self, fname, mode='r', **kwargs):
        Dataset.__init__(self, filename=fname, mode=mode, **kwargs)
#        print "Opening file %s" % fname
#        print "Done"

#        print "Initing map"
        self.proj = Proj(
                proj='lcc',
                lat_0 = self.YCENT,
                lon_0 = self.XCENT,
                lat_1 = self.P_ALP, 
                lat_2 = self.P_BET,
                a=6370000.,
                b=6370000.
                )
#        print "Done"

    def __setattr__(self, name, value):
        try:
            Dataset.__setattr__(self, name, value)
        except TypeError:
            self.__dict__[name] = value
    
    def lonlat_to_xy(self, lon, lat):
        return self.proj(lon,lat)

    def lonlat_to_ij(self, lon, lat, dot=False):
        x, y = self.lonlat_to_xy(lon, lat)
        i = (x - self.XORIG) / self.XCELL
        j = (y - self.YORIG) / self.YCELL
        if not dot:
            i -= 0.5
            j -= 0.5
        return i, j
    
    def xy_to_lonlat(self, x, y):
        return self.proj(x, y, inverse=True)

    def ij_to_xy(self, i, j, dot=False):
        if not dot:
            i = i + 0.5
            j = j + 0.5
        x = i * self.XCELL + self.XORIG
        y = j * self.YCELL + self.YORIG
        return x, y

    def ij_to_lonlat(self, i, j, **kwargs):
        x, y = self.ij_to_xy(i, j, **kwargs)
        return self.xy_to_lonlat(x, y)

#    def get_data_pos(self, lon, lat, height, dot=False):
#        i, j = self.lonlat_to_ij(lon, lat, dot=dot)
#        if i < -0.5 or j < -0.5 \
#                or i > self.NCOLS - 0.5 \
#                or j > self.NROWS - 0.5:
#            return None
#        elif i < 0.5 or j < 0.5 \
#                or i > self.NCOLS - 1.5 \
#                or j > self.NROWS - 1.5:
#            pass       
#
#    def get_data(self, varname, data_pos):
#        pass

    def get_time_serie(self, varname, lon, lat, level=0, **kwargs):
        i, j = self.lonlat_to_ij(lon, lat, **kwargs)
        i = int(round(i))
        j = int(round(j))
        return self.variables[varname][:,level, j, i]
    
    def get_xy_arrays(self, **kwargs):
        i = np.arange(self.NCOLS, dtype=int)
        j = np.arange(self.NROWS, dtype=int)
        ii, jj = np.meshgrid(i, j)
        x, y = self.ij_to_xy(ii, jj, **kwargs)
        return x, y
    
    def get_lonlat_arrays(self, **kwargs):
        x, y = self.get_xy_arrays(**kwargs)
        lon, lat = self.xy_to_lonlat(x, y)
        return lon, lat
    
