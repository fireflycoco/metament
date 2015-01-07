#!/usr/bin/env python

# lidar.py
"""class LidarDataset is a container for lidar data.
    TODO: add support for loading short lines.
    TODO: make 'must'/'optional' work
"""
import os, sys
from datetime import datetime, timedelta
from copy import copy, deepcopy
import numpy as np
from netCDF4 import Dataset, date2num, num2date
from metlib.datetime.datetime_bin import datetime_bin

__all__ = ['LidarDataset']
_std_datetime_fmt = "%Y-%m-%d %H:%M:%S"
_std_datetime_units = "seconds since 1970-01-01 00:00:00"
_NO_DESC_STR = "No description available"
# # Table driven:
# # varname, format, dim shape, 'must'/'optional', average method

_special_avermethods = {
        'data':'sum',
        'datetime':'special',
        'shots_sum':'sum',
        'dead_time_corrected':'first',
        'background':'sum',
        'energy':'sum',
        'background_std_dev':'sqr_mean_sqrt'
        }
_important_attrs = set([
    'lidarname', 'bin_time', 'bin_size', 
    'first_data_bin', 'start_datetime', 'end_datetime', 
    'number_records', 'number_bins', 'number_channels',
    'max_range',
    'data_aver_method','desc'
    ])

_important_vars = set([
    'data', 'energy', 'background',
    'datetime','distance',
    ])

class LidarDataset(object):
    """Representing Lidar's Dataset"""
    def __init__(self, fnames, bin_num=None, **kwargs):
        """doc
        kwargs:
            fnames: a seq of filenames or a single filename or a str of filenames seperated with comma.
            bin_num: clipping in BIN dimesion when loading files.
        """
        self.init_clean()
        if isinstance(fnames, np.string_):
            fnames = str(fnames)
        if isinstance(fnames, (str, unicode)):
            fnames = fnames.split(',')
        try:
            self.read_one_file(fnames[0], **kwargs)
        except RuntimeError:
            if len(fnames) == 1:
                raise
            else:
                self.__init__(fnames[1:], bin_num=bin_num, **kwargs)
                return
        if len(fnames) > 1:
            self.append_files(fnames[1:], **kwargs)

        if bin_num is not None:
            self.resize_bin(0,bin_num)

    def append_datasets(self, datasets):
        """append one or more LidarDatasets to this dataset.
        datasets: a seq of datasets or a single LidarDataset object"""
        if len(datasets) == 0:
            return
        if isinstance(datasets, LidarDataset):
            datasets = [datasets]
        for vname in self.vars:
            dims = self.var_dims[vname]
            if 'TIME' in dims:
                stack_func = np.hstack if len(dims) == 1 else np.vstack
                self.vars[vname] = stack_func(tuple([self[vname]] + [d[vname] for d in datasets]))
        self.recheck_time()
        del datasets

    def init_clean(self):
        """clean up everything"""
        self.dims = {}
        self.vars = {}
        self.attrs = {}
        self.var_dims = {}
        self.var_aver_methods = {}

    def read_one_file(self, fname, **kwargs):
        """read one file for initing self"""
        self.init_clean()
        f = Dataset(fname, **kwargs)
        for d in f.dimensions:
            self.dims[str(d)] = len(f.dimensions[d])
        if 'TIME' not in self.dims or self.dims['TIME'] == 0:
            raise RuntimeError("%s contains no data" % fname)
        for v in f.variables:
            strv = str(v)
            ncv = f.variables[v]
            v_arr = ncv[:]
            # datetime handling
            if v == 'datetime':
                datetime_type = type(v_arr[0])
                if datetime_type in (np.string_ , str, unicode):
                    v_arr = np.array([datetime.strptime(datestr, _std_datetime_fmt) for datestr in v_arr])
                else:
                    v_arr = num2date(v_arr, units=_std_datetime_units)
            self.vars[strv] = v_arr
            self.var_dims[strv] = tuple([str(dimname) for dimname in ncv.dimensions])
            # aver_method
            if 'aver_method' in ncv.ncattrs():
                self.var_aver_methods[strv] = str(ncv.getncattr('aver_method'))
            elif v in _special_avermethods:
                self.var_aver_methods[strv] = _special_avermethods[v]
            else:
                self.var_aver_methods[strv] = 'mean'
        for a in f.ncattrs():
            self.attrs[str(a)] = f.getncattr(a)
        self.recheck_time()
        f.close()
    
    def append_files(self, fnames):
        """append one or more files to the dataset
        fnames: a seq of filenames or a single filename or a str of filenames seperated with comma."""
        if isinstance(fnames, np.string_):
            fnames = str(fnames)
        if isinstance(fnames, (str, unicode)):
            fnames = fnames.split(',')
        tmpd = []
        for fn in fnames:
            try:
                d = LidarDataset(fn)
                tmpd.append(d)
            except RuntimeError:
                pass
        self.append_datasets(tmpd)
        del tmpd

    def save(self, fname, use_datetime_str=True):
        """Save into a netCDF4 file"""
        if os.path.exists(fname):
            if os.path.islink(fname):
                os.unlink(fname)
            else:
                os.remove(fname)
        f = Dataset(fname, 'w', format='NETCDF4')
        for dname, length in self.dims.items():
            if dname == 'TIME':
                f.createDimension(dname)
            else:
                f.createDimension(dname, length)
        for vname in self.vars:
            dimnames = self.var_dims[vname]
            t = self.vars[vname].dtype.str.lstrip('<>|=')
            if type(vname) is unicode:
                vname = str(vname)
            if vname == 'datetime':
                if use_datetime_str:
                    f.createVariable(vname, str, dimnames)
                    for i, dt in enumerate(self.vars[vname]):
                        f.variables[vname][i] = self.vars[vname][i].strftime(_std_datetime_fmt)
                else:
                    f.createVariable(vname, 'i4', dimnames)
                    f.variables[vname].setncattr('units', _std_datetime_units)
                    f.variables[vname][:] = date2num(self.vars[vname], _std_datetime_units)
            else:
                if t.startswith('O'):
                    t = str
                f.createVariable(vname, t, dimnames)
                f.variables[vname][:] = self.vars[vname]
            f.variables[vname].setncattr('aver_method', self.var_aver_methods[vname])
    
        for attr in self.attrs:
            if type(attr) is unicode:
                attr = str(attr)
            if isinstance(self.attrs[attr], datetime):
                f.setncattr(attr, self.attrs[attr].strftime(_std_datetime_fmt))
            else:
                a = str(self.attrs[attr]) if type(self.attrs[attr]) is unicode else self.attrs[attr]
                f.setncattr(attr, a)
        f.close()

    def time_average(self, tdelta, starttime=None, endtime=None):
        """Averaging in TIME dimesion.
        tdelta is a datetime.timedelta object or an int of minutes.
        starttime and endtime works in the same way as metlib.datetime.datetime_bin.datetime_bin .

        Notice: Some variables are summed instead of averaged. Including:
        data, energy, background, etc.
        """
        if type(tdelta) is not timedelta:
            tdelta = timedelta(minutes=tdelta)
        dts = self.vars['datetime']
        tbins, info = datetime_bin(dts, tdelta, starttime=starttime, endtime=endtime, return_bin_info=True)
        new_time_steps = len(tbins)
        tmp = dict()

        for vname in self.vars:
            dimnames = self.var_dims[vname]
            if 'TIME' in dimnames:
                new_shape = list(self.vars[vname].shape)
                new_shape[0] = new_time_steps
                tmp[vname] = np.zeros(tuple(new_shape), dtype=self.vars[vname].dtype)

        for i, w in enumerate(tbins):
            for vname in self.vars:
                dimnames = self.var_dims[vname]
                aver_method = self.var_aver_methods[vname]
                t = self.vars[vname].dtype.char
                if 'TIME' not in dimnames:
                    continue
                if vname == 'datetime':
                    # # using each tbins' starttime as datetime
                    tmp[vname][i] = info[i][0]
                elif len(w[0]) == 0:
                    # # filling with nan or 0
                    if t not in ('f', 'd', 'g'):
                        tmp[vname][i] = 0
                    else:
                        tmp[vname][i] = np.nan
                else:
                    if aver_method == 'sum':
                        tmp[vname][i] = np.nansum(self.vars[vname][w], axis=0)
                    elif aver_method == 'mean':
                        tmp[vname][i] = np.ma.masked_invalid(self.vars[vname][w]).mean(axis=0)
                    elif aver_method == 'first':
                        tmp[vname][i] = self.vars[vname][w][0]
                    elif aver_method == 'sqr_mean_sqrt':
                        # # for std_dev 
                        tmp[vname][i] = np.sqrt(np.mean(np.ma.masked_invalid(self.vars[vname][w]) ** 2, axis=0))
                    elif aver_method == 'min':
                        seled = self.vars[vname][w]
                        maskarr = np.ma.masked_invalid(seled)
                        if maskarr.count() == 0:
                            the_min = seled[0]
                        else:
                            the_min = maskarr.min(axis=0)
                    elif aver_method == 'max':
                        seled = self.vars[vname][w]
                        maskarr = np.ma.masked_invalid(seled)
                        if maskarr.count() == 0:
                            the_max = seled[0]
                        else:
                            the_max = maskarr.max(axis=0)
                        tmp[vname][i] = the_max
                    elif aver_method == 'positive_min':
                        seled = self.vars[vname][w]
                        maskarr = np.ma.masked_where(-(seled >= 0), seled)
                        if maskarr.count() == 0:
                            the_min = seled[0]
                        else:
                            the_min = maskarr.min(axis=0)
                        tmp[vname][i] = the_min
                    elif aver_method == 'positive_max':
                        seled = self.vars[vname][w]
                        maskarr = np.ma.masked_where(-(seled <= 0), seled)
                        if maskarr.count() == 0:
                            the_max = seled[0]
                        else:
                            the_max = maskarr.max(axis=0)
                        tmp[vname][i] = the_max

        self.vars.update(tmp)
        self.recheck_time()

    def __len__(self):
#        return self.dims['TIME']
        return self.vars['data'].shape[0]

    def recheck_time(self):
        """recheck stuffs about time dimension"""
        n_tbins = len(self.vars['datetime'])
        self.attrs['start_datetime'] = self.vars['datetime'][0]
        if len(self) >= 2:
            self.attrs['end_datetime'] = self.vars['datetime'][-1] + (self.vars['datetime'][-1] - self.vars['datetime'][-2])
        elif self.vars['trigger_frequency'][-1] != 0:
            self.attrs['end_datetime'] = self.vars['datetime'][-1] + timedelta(seconds = int(self.vars['shots_sum'][-1] / self.vars['trigger_frequency'][-1]))
        self.attrs['number_records'] = n_tbins
        self.dims['TIME'] = n_tbins

    def get_timedeltas(self, return_seconds=True):
        """Get array of each record's time span. 
        if return_seconds is True: return in seconds
        else: return in timedelta objects.
        """
        seconds = np.zeros(self.dims['TIME'], dtype='i8')
        w = np.where(self.vars['trigger_frequency'] > 0)
        seconds[w] = self.vars['shots_sum'][w] / self.vars['trigger_frequency'][w]
        if return_seconds:
            return seconds
        else:
            tds = np.array([timedelta(seconds = int(s)) for s in seconds])
            return tds

    def get_end_datetimes(self):
        """Get array of each record's end datetimes"""
        tds = self.get_timedeltas(return_seconds=False)
        return self.vars['datetime'] + tds

    def get_mid_datetimes(self):
        """Get array of each record's mid datetimes"""
        tds = self.get_timedeltas(return_seconds=False)
        return self.vars['datetime'] + tds / 2

    def keep_indice(self, indice):
        """Keep only data in the indice (Time dimension)"""
        # FIXME when only one record is kept, there's a bug
        where_to_keep = (np.array(indice),)
        tmp = dict()
        for vname in self.vars:
            if 'TIME' in self.var_dims[vname]:
                tmp[vname] = self.vars[vname][where_to_keep]
        self.vars.update(tmp)
        self.recheck_time()

    def drop_indice(self, indice):
        """Drop unwanted data in the indice(Time dimension)"""
        my_indice = set(range(len(self)))
        my_indice.difference_update(set(indice))
        to_keep = list(my_indice)
        to_keep.sort()
        self.keep_indice(to_keep)
    
    
    def drop_periods(self, periods):
        """Drop unwanted data in the time periods.
        periods:
            A list of tuples of (beg_time, end_time)
        """
        to_drop = []
        dts = self.vars['datetime']
        for (beg, end) in periods:
            w = np.where((dts>=beg) & (dts<end))
            to_drop.extend(list(w[0]))
        self.drop_indice(to_drop)

    def keep_periods(self, periods):
        """Keep the periods only, drop other periods.
        periods:
            A list of tuples of (beg_time, end_time)
        """
        to_keep = []
        dts = self.vars['datetime']
        for (beg, end) in periods:
            w = np.where((dts>=beg) & (dts<end))
            to_keep.extend(list(w[0]))
        self.keep_indice(to_keep)

    def trim_period(self, starttime, endtime):
        """Trim unwanted data, leaving data between start_time and end_time only.
        """
        self.keep_periods([(starttime, endtime)])

    def __getitem__(self, key):
        """return a new LidarDataset object that contains the slice (in TIME dimension) if key is slice or integer. 
        return corresponding data if key is str
        """
        if isinstance(key, slice):
            new_data = self.copy()
            for vname in self.vars:
                if 'TIME' not in self.var_dims[vname]:
                    continue
                new_data.vars[vname] = new_data.vars[vname][key]
            new_data.recheck_time()
            return new_data
        elif isinstance(key, (int, long, np.integer)):
            new_data = self.copy()
            for vname in self.vars:
                if 'TIME' not in self.var_dims[vname]:
                    continue
                new_data.vars[vname] = np.array(new_data.vars[vname][key])[np.newaxis, ...]
            new_data.recheck_time()
            return new_data
        elif isinstance(key, (str, unicode)):
            if key in self.vars:
                return self.vars[key]
            elif key in self.attrs:
                return self.attrs[key]

    def __setitem__(self, key, value):
        if isinstance(key, (str, unicode)):
            if key in self.vars:
                self.vars[key][:] = value
            else:
                self.attrs[key] = value
        else:
            sys.stderr.write('Warning: key is not str, not setting LidarDataset property\n')

    def __delitem__(self, key):
        if isinstance(key, (str, unicode)):
            if key in self.attrs:
                if key in _important_attrs:
                    sys.stderr.write("Warning: Not deleting important attr: %s \n" % key)
                    return
                try:
                    del self.attrs[key]
                except:
                    sys.stderr.write("Warning: cannot delete attr: %s \n" % key)
            elif key in self.vars:
                self.del_var(key)
            else:
                sys.stderr.write("Warning: no such attr: %s \n" % key)

    def __contains__(self, key):
        return key in self.attrs or key in self.vars
    
    def add_var(self, vname, dims, dtype='f4', aver_method='mean'):
        """Add new var into self.vars"""
        if vname in self.vars:
            sys.stderr.write("Warning: %s already in the dataset. Not adding\n" % vname)
            return
        # TODO: check dimensions
        self.vars[vname] = np.zeros(tuple([self.dims[d] for d in dims]), dtype=dtype)
        self.var_dims[vname] = dims
        self.var_aver_methods[vname] = aver_method

    def del_var(self, vname):
        """Del var in self.vars.
        """
        if vname not in self.vars:
            sys.stderr.write("Warning: %s not in the dataset\n" % vname)
            return
        if vname in _important_vars:
            sys.stderr.write("Warning: %s is important, thus not deleting\n" % vname)
            return
        del self.vars[vname]
        del self.var_dims[vname]
        del self.var_aver_methods[vname]

    def resize_bin(self, start_i, end_i):
        """resize data's BIN dim in python's manner: data[start_i:end_i]
        """
        if end_i > self.attrs['number_bins']:
            pass
        else:
            self.vars['data'] = self.vars['data'][:,:,start_i:end_i]
            self.vars['distance'] = self.vars['distance'][start_i:end_i]
            self.attrs['number_bins'] = np.max(end_i - start_i, 0)
            self.attrs['first_data_bin'] -= start_i

            self.dims['BIN'] = self.attrs['number_bins']
    
    def copy(self):
        return deepcopy(self)
   
    def __add__(self, other):
        res = self.copy()
        res.append_datasets(other)
        return res

    def __str__(self):
        if 'desc' not in self.attrs:
            desc = _NO_DESC_STR
        else:
            desc = self.attrs['desc']
        return """    lidarname: %s
    desc: %s
    time period: %s - %s
    dims: %s
    vars: %s
    attrs: %s
    """ % ( self.attrs['lidarname'],
            desc,
            self.attrs['start_datetime'], self.attrs['end_datetime'],
            self.dims, 
            self.vars.keys(),
            self.attrs.keys()
            )

    def __repr__(self):
        return self.__str__()
    
    @property
    def desc(self):
        return self.attrs['desc']
    @desc.setter
    def desc(self, value):
        self.attrs['desc'] = value

if __name__ == '__main__':
    d = LidarDataset('00.nc,01.nc,02.nc,03.nc')
    print d
