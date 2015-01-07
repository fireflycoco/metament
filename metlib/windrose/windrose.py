#!/usr/bin/env python
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from mpl_toolkits.basemap import Basemap

__all__ = ['WindroseBase', 'Windrose']
# <codecell>

_WindroseDefaults = dict(
     hollow=0.15,
     rmin=0.0,
     rmax=10.0,
     rstep=2.0,
     rlabeldir=60.0,
     rlabelfmt='%.1f',
     cmap=plt.cm.jet,
     gridcolor='black',
     gridlinewidth=1.0,
     griddashes=[1,1],
     circlecolor='black',
     circlelinewidth=1.0,
     bgcolor=None,
     griddirnum=8,
     rlabelsize='small',
     rlabelweight='normal',
     rlabelcolor='black',
     dirlabelsize='medium',
     dirlabelweight='semibold',
     dirlabelcolor='black',
     )

class WindroseBase(object):
    def __init__(self, wd, ws, data, ax=None, **kwargs):
        self.__dict__.update(_WindroseDefaults)
        self.wd = np.copy(wd)
        self.ws = np.copy(ws)
        self.data = np.copy(data)
        self.ax = plt.gca() if ax is None else ax
        self.fig = self.ax.figure
        
        self.m = Basemap(projection='ortho', 
            lat_0=-90, lon_0=0, 
            resolution='c',
            )
        
        for k, v in kwargs.iteritems():
            if k in _WindroseDefaults:
                self.__dict__[k] = v
            else:
                sys.stderr.write('Warning: kwarg: %s does not fit for Windrose\n' % k)
        self.set_rlim()
        self.set_griddirnum()

    def set_rlim(self, rmin=None, rmax=None, rstep=None, hollow=None):
        if rmin is not None:
            self.rmin = rmin
        if rmax is not None:
            self.rmax = rmax
        if rstep is not None:
            self.rstep = rstep
        if hollow is not None:
            self.hollow = hollow
        self._grid_rs = np.arange(self.rmin, self.rmax+self.rstep/10.0, self.rstep)
        self._par_lats = self._r2lat(self._grid_rs)
        self._latmax = np.max(np.abs(self._par_lats))
    
    def set_griddirnum(self, griddirnum=None):
        if griddirnum is not None: 
            self.griddirnum=griddirnum
        self._mer_lons= np.arange(0, 360.0, 360.0/self.griddirnum)
    
    def draw(self):
        plt.sca(self.ax)
        self.m.drawparallels(self._par_lats, latmax=self._latmax, 
            linewidth=self.gridlinewidth, color=self.gridcolor, dashes=self.griddashes)
        self.m.drawmeridians(self._mer_lons, latmax=self._latmax,
            linewidth=self.gridlinewidth, color=self.gridcolor, dashes=self.griddashes)
        self.m.drawmapboundary(color=self.circlecolor, linewidth=self.circlelinewidth, fill_color=self.bgcolor)
        self.ax.text(0.5, 1.03, 'N', ha='center', va='bottom', transform=self.ax.transAxes, fontsize=self.dirlabelsize, fontweight=self.dirlabelweight, color=self.dirlabelcolor)
        self.ax.text(0.5, -0.03, 'S', ha='center', va='top', transform=self.ax.transAxes, fontsize=self.dirlabelsize, fontweight=self.dirlabelweight, color=self.dirlabelcolor)
        self.ax.text(1.03, 0.5, 'E', ha='left', va='center', transform=self.ax.transAxes, fontsize=self.dirlabelsize, fontweight=self.dirlabelweight, color=self.dirlabelcolor)
        self.ax.text(-0.03, 0.5, 'W', ha='right', va='center', transform=self.ax.transAxes, fontsize=self.dirlabelsize, fontweight=self.dirlabelweight, color=self.dirlabelcolor)
        for r, pl in zip(self._grid_rs, self._par_lats):
            lon, lat = self.m(self.rlabeldir, pl)
            self.ax.text(lon, lat, self.rlabelfmt % r, 
                color=self.rlabelcolor,
                fontsize=self.rlabelsize, 
                )
                

    def _r2lat(self, r):
        bigR = (self.rmax - self.rmin) / (1.0 - self.hollow)
        real_r = (bigR*self.hollow + np.array(r)) / bigR
        lat = -np.rad2deg(np.arccos(real_r))
        return lat
    
    def extend_wd(self):
        return np.hstack((self.wd, self.wd[-1] * 2 -self.wd[-2]))
    
    def extend_ws(self):
        return np.hstack((self.ws, self.ws[-1] * 2 - self.ws[-2]))
    
    def extend_data_for_contour(self):
        return np.hstack((self.data, self.data[:, 0:1]))
    
    def offset_data(self):
        wd_step = self.wd[1] - self.wd[0]
        self.wd -= wd_step/2

# <codecell>

class Windrose(WindroseBase):
    def __init__(self, wd, ws, data, ax=None, offset_data=False, **kwargs):
        WindroseBase.__init__(self, wd, ws, data, ax=ax, **kwargs)
        self.ext_wd = self.extend_wd()
        self.ext_ws = self.extend_ws()
        if offset_data:
            self.offset_data()
        
        
    def pcolor(self, **kwargs):
        WindroseBase.draw(self)
        NT, NR = np.meshgrid(self.ext_wd, self.ext_ws)
        lon, lat = self.m(NT, self._r2lat(NR))
        im = self.m.pcolormesh(lon, lat, self.data, **kwargs)
        return im

    def contourf(self, **kwargs):
        ext_data = self.extend_data_for_contour()
        WindroseBase.draw(self)
        NT, NR = np.meshgrid(self.ext_wd, self.ws)
        lon, lat = self.m(NT, self._r2lat(NR))
        cs = self.m.contourf(lon, lat, ext_data, **kwargs)
        return cs
    
    def contour(self, **kwargs):
        ext_data = self.extend_data_for_contour()
        WindroseBase.draw(self)
        NT, NR = np.meshgrid(self.ext_wd, self.ws)
        lon, lat = self.m(NT, self._r2lat(NR))
        cs = self.m.contour(lon, lat, ext_data, **kwargs)
        return cs
    

# <codecell>
if __name__ == '__main__':
    delta = 0.025
    delta = 0.25
    x = np.arange(-3.0, 3.0, delta)
    y = np.arange(-2.0, 2.0, delta)
    X, Y = np.meshgrid(x, y)
    Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    # difference of Gaussians
    Z = 10.0 * (Z2 - Z1)

    T = (x+3.0) * 60
    R = (y+2.0) * 2.5

    NT = (X+3.0) * 60
    NR = (Y+2.0) * 2.5

    # <codecell>

    wb = WindroseBase(T, R, Z, bgcolor='c')
    wb.draw()
    plt.show()

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_axes((0.1, 0.1, 0.7, 0.8))
    cax = fig.add_axes((0.85, 0.1, 0.1, 0.8))
    wr = Windrose(T, R, Z, ax=ax, hollow=0.1, circlecolor='r', circlelinewidth=3, offset_data=True)
#    im = wr.contourf(cmap=plt.cm.Pastel2)
    im = wr.pcolor(cmap=plt.cm.Pastel2)
    plt.colorbar(im, cax)
    plt.show()
    wr.fig.savefig('./test_windrose.png')


