#!/usr/bin/env python2.7

# wps.py
# FIXME: Debugging!

import os, sys
from StringIO import StringIO
import numpy as np
import pyproj
from metlib.data.grid import XYGrid2D

#from netCDF4 import Dataset

def parse_wps_namelist(namelist):
    """parse_wps_namelist() parses wps namelist.
    namelist:
        a namelist file or a str of (part of) namelist.

    Returns:
        a dict that contains namelist info.
    """

    if '\r' in namelist or '\n' in namelist:
        f = StringIO(namelist)
    else:
        f = open(namelist)

    result = {'default':dict() }
    current_group = 'default'
    for line in f:
        line = line.strip()
        if len(line) == 0 or line == '/':
            continue
        if line.startswith('&'):
            current_group = line[1:]
            if current_group not in result:
                result[current_group] = dict()
            continue
#        print line
        varname, values = line.split('=')
        varname = varname.strip()
        values = values.strip()
        value_list = values.split(',')
        value_list = [value.strip() for value in value_list]
        converted_value_list = []
        for value in value_list:
            if value == "":
                continue
            if (value[0] == '"' and value[-1] == '"') or \
                    (value[0] == "'" and value[-1] == "'"):
                converted_value_list.append(value[1:-1])
                continue
            try:
                converted_value = int(value)
            except ValueError:
                try:
                    converted_value = float(value)
                except ValueError:
                    converted_value = value
            finally: 
                converted_value_list.append(converted_value)
        result[current_group][varname] = converted_value_list

    return result

class Domainer(object):
    def __init__(self, namelist):
        if not isinstance(namelist, dict):
            namelist = parse_wps_namelist(namelist)
        self.geogrid = namelist['geogrid']
        self.domains = range(1, len(self.geogrid['parent_id'])+1)
        self.projections = dict()
        self.truelat1 = self.geogrid['truelat1'][0]
        self.truelat2 = self.geogrid['truelat2'][0]
        self.ref_lat = self.geogrid['ref_lat'][0]
        self.ref_lon = self.geogrid['ref_lon'][0]
        self.stand_lon = self.geogrid['stand_lon'][0]
        
        if self.stand_lon != self.ref_lon:
            raise ValueError("Standard Longitude must equal Ref longitude!")

        self.dx = float(self.geogrid['dx'][0])
        self.dy = float(self.geogrid['dy'][0])
        self.earth_r_a = 6370000.0
        self.earth_r_b = 6370000.0
        self.proj = pyproj.Proj(
                proj = 'lcc',
                lat_0 = self.ref_lat,
                lon_0 = self.ref_lon,
                lat_1 = self.truelat1,
                lat_2 = self.truelat2,
                a = self.earth_r_a,
                b = self.earth_r_b,
                )

        self.domain_dx = dict()
        self.domain_dy = dict()
        self.xorig = dict()
        self.yorig = dict()
        self.xgridnum = dict()
        self.ygridnum = dict()

        for i in self.domains:
            self.xgridnum[i] = self.geogrid['e_we'][i-1]
            self.ygridnum[i] = self.geogrid['e_sn'][i-1]

        for i in self.domains:
            family_tree = [i]
            parent = self.geogrid['parent_id'][family_tree[-1]-1]
            while family_tree[-1] != 1:
                family_tree.append(parent)
                parent = self.geogrid['parent_id'][family_tree[-1]-1]
            family_tree.reverse()

            self.xorig[i] = 0.0
            self.yorig[i] = 0.0
            for member in family_tree:
                if member == 1:
                    cur_dx = self.dx
                    cur_dy = self.dy
                    self.xorig[i] -= ((self.xgridnum[member]-1) / 2.0 ) * cur_dx
                    self.yorig[i] -= ((self.ygridnum[member]-1) / 2.0 ) * cur_dy
                else:
                    self.xorig[i] += float(self.geogrid['i_parent_start'][member-1] - 1) * cur_dx 
                    self.yorig[i] += float(self.geogrid['j_parent_start'][member-1] - 1) * cur_dy 
                cur_dx /= float(self.geogrid['parent_grid_ratio'][member-1])
                cur_dy /= float(self.geogrid['parent_grid_ratio'][member-1])

            self.domain_dx[i] = cur_dx
            self.domain_dy[i] = cur_dy

        self.dot_grids = dict()
        self.cross_grids = dict()
        for i in self.domains:
            dot_IX, dot_JY = np.meshgrid(np.arange(self.xgridnum[i]), np.arange(self.ygridnum[i]))
            cross_IX, cross_JY = np.meshgrid(np.arange(self.xgridnum[i]-1), np.arange(self.ygridnum[i]-1))
            dot_lons, dot_lats = self.ij_to_lonlat(dot_IX, dot_JY, domain=i, dot=True)
            cross_lons, cross_lats = self.ij_to_lonlat(cross_IX, cross_JY, domain=i, dot=False)
            self.dot_grids[i] = XYGrid2D(dot_lons, dot_lats)
            self.cross_grids[i] = XYGrid2D(cross_lons, cross_lats)

    def lonlat_to_xy(self, lon, lat):
        return self.proj(lon, lat)

    def lonlat_to_ij(self, lon, lat, domain=1, dot=False):
        x, y = self.lonlat_to_xy(lon, lat)
        i = (x - self.xorig[domain]) / self.domain_dx[domain]
        j = (y - self.yorig[domain]) / self.domain_dy[domain]
        if not dot:
            i = i - 0.5
            j = j - 0.5
        return i, j
    
    def xy_to_lonlat(self, x, y):
        return self.proj(x, y, inverse=True)

    def ij_to_xy(self, i, j, domain=1, dot=False):
        if not dot:
            i = i + 0.5
            j = j + 0.5
        x = i * self.domain_dx[domain] + self.xorig[domain]
        y = j * self.domain_dy[domain] + self.yorig[domain]
        return x, y

    def ij_to_lonlat(self, i, j, domain=1, dot=False):
        x, y = self.ij_to_xy(i, j, domain=domain, dot=dot)
        return self.xy_to_lonlat(x, y)

if __name__ == '__main__':
    test_nl = """
&share
 wrf_core = 'ARW',
 max_dom = 2,
 start_date = '2008-03-24_12:00:00','2008-03-24_12:00:00',
 end_date   = '2008-03-24_18:00:00','2008-03-24_12:00:00',
 interval_seconds = 21600,
 io_form_geogrid = 2
/

&geogrid
 parent_id         =   1,   1,  2,
 parent_grid_ratio =   1,   3,  3,
 i_parent_start    =   1,  31,  20,
 j_parent_start    =   1,  17,  15,
 s_we              =   1,   1, 1
 e_we              =  74, 112, 100,
 s_sn              =   1,   1, 1,
 e_sn              =  61,  97, 80,
 geog_data_res     = '10m','2m','30s'
 dx = 30000,
 dy = 30000,
 map_proj = 'lambert',
 ref_lat   = 34.83,
 ref_lon   = -81.03,
 truelat1  =  30.0,
 truelat2  =  60.0,
 stand_lon = -98.,
 geog_data_path = '/mmm/users/wrfhelp/WPS_GEOG/'
/
"""

    res = parse_wps_namelist(test_nl)

    print res

    dmr = Domainer(res)

