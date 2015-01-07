#!/usr/bin/env python

# misc.py

import os
import re
#from datetime import datetime, timedelta
#from dateutil.parser import parse
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
#from matplotlib import mlab
#from netCDF4 import Dataset

__all__ = ['axgrid', 'axstack']

def axgrid(fig=None, nrow=2, ncol=2, left=0.1, right=0.1, top=0.1, bottom=0.1, hspace=0.05, vspace=0.05, no_extra_xticklabels=True, no_extra_yticklabels=True, sharex=False, sharey=False, sharexy=False, **kwargs):
    """axgrid creates a grid of axes on a fig.
    Parameters:
        fig: the fig. if fig is None: fig = plt.gcf()
        nrow: number of rows.
        ncol: number of columns.
        left: left edge width (0-1).
        right: right edge width (0-1).
        top: top edge width (0-1).
        bottom: bottom edge width (0-1).
        hspace: horizontal space between axes.
        vspace: vertical space between axes.
        no_extra_xticklabels: if True, only xticklabels on the bottom row are presented.
        no_extra_yticklabels: if True, only yticklabels on the left row are presented.
        sharex: if True, axes on each column will share the same xaxis.
        sharey: if True, axes on each row will share the same yaxis.
        sharexy: if True, all axes will share the same xaxis and yaxis.
        kwargs: pass to fig.add_axes()
    Returns:
        an np.ndarray of axes. 
    """
    if fig is None:
        fig = plt.gcf()
    thegrid = np.zeros((nrow, ncol), dtype='O')

    w = (1.0 - left - right - hspace*(ncol-1)) / ncol
    h = (1.0 - top - bottom - vspace*(nrow-1)) / nrow

    for idx in np.ndindex(nrow, ncol):
        pos = (left+idx[1]*(w+hspace), 1.0-top-h-idx[0]*(h+vspace), w, h)
        sharexyd = dict()
        if sharexy:
            if idx[0] >= 1 or idx[1] >= 1:
                sharexyd['sharex'] = thegrid[0, 0]
                sharexyd['sharey'] = thegrid[0, 0]
        else:
            if sharex and idx[0] >= 1:
                sharexyd['sharex'] = thegrid[0, idx[1]]
            if sharey and idx[1] >= 1:
                sharexyd['sharey'] = thegrid[idx[0], 0]
        final_kws = sharexyd.copy()
        final_kws.update(kwargs)
        thegrid[idx] = fig.add_axes(pos, **final_kws)

    if no_extra_xticklabels:
        for ax in thegrid[:-1, :].flat:
            plt.setp(ax.get_xticklabels(), visible=False)
    if no_extra_yticklabels:
        for ax in thegrid[:, 1:].flat:
            plt.setp(ax.get_yticklabels(), visible=False)
    return thegrid

def axstack(fig=None, pos=(0.1, 0.1, 0.8, 0.8), n=3, xoffset=0.1, yoffset=0.1, axis_color='k', axis_linewidth=1.0, tick_len=2.0, no_extra_xticklabels=True, no_extra_yticklabels=True, sharexy=True, perspective=0.9, **kwargs):
    def map_pos(x_in, y_in):
        return x_in*pos[2] + pos[0], y_in*pos[3] + pos[1]

    if fig is None:
        fig = plt.gcf()
    axs = np.zeros(n+1, dtype='O')

    axframe = fig.add_axes(pos)
    plt.setp(axframe.get_xticklines(), visible=False)
    plt.setp(axframe.get_yticklines(), visible=False)
    plt.setp(axframe.get_xticklabels(), visible=False)
    plt.setp(axframe.get_yticklabels(), visible=False)
    for spine in axframe.spines.values():
        spine.set_visible(False)

    axs[n] = axframe

    if np.isscalar(xoffset):
        xoffset = np.hstack((0.0, np.add.accumulate([xoffset] * (n-1))))
    if np.isscalar(yoffset):
        yoffset = np.hstack((0.0, np.add.accumulate([yoffset] * (n-1))))

    if np.isscalar(axis_color):
        axis_color = [axis_color] * n

    WW, HH = pos[2:4]
    wws = np.zeros(n)
    hhs = np.zeros(n)
    wws[-1], hhs[-1] = 1.0-xoffset[-1], 1.0-yoffset[-1]
    pers_ratios = np.interp(xoffset, [xoffset[0], xoffset[-1]], [1.0/ perspective, 1.0])
    inv_pers_ratios = pers_ratios / pers_ratios[0]
    wws[:] = wws[-1] * pers_ratios
    hhs[:] = hhs[-1] * pers_ratios
    if np.isscalar(axis_linewidth):
        axis_linewidth = axis_linewidth * inv_pers_ratios
    if np.isscalar(tick_len):
        tick_len = tick_len * inv_pers_ratios

    for i in range(n-1, -1, -1):
        x0, y0 = map_pos(xoffset[i], yoffset[i])
        if sharexy:
            ax = fig.add_axes((x0, y0, wws[i]*WW, hhs[i]*HH), axis_bgcolor='none', sharex=axframe, sharey=axframe, **kwargs)
        else:
            ax = fig.add_axes((x0, y0, wws[i]*WW, hhs[i]*HH), axis_bgcolor='none', **kwargs)
        plt.setp(ax.get_xticklines(), color=axis_color[i])
        plt.setp(ax.get_yticklines(), color=axis_color[i])
        ax.tick_params(length=tick_len[i]) 
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        for direction in ['left', 'bottom']:
            ax.spines[direction].set_linewidth(axis_linewidth[i])
            ax.spines[direction].set_color(axis_color[i])
        for direction in ['right', 'top']:
            ax.spines[direction].set_visible(False)
        axs[i] = ax
    
    for ax in axs[1:n]:
        if no_extra_xticklabels:
            plt.setp(ax.get_xticklabels(), visible=False)
        if no_extra_yticklabels:
            plt.setp(ax.get_yticklabels(), visible=False)

    return axs

if __name__ == '__main__':
    pass
