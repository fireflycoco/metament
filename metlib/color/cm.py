#!/usr/bin/env python

from matplotlib.colors import LinearSegmentedColormap as LSCM
from metlib.color.color_collections import *

_bb = 'CornflowerBlue'
_br = 'DarkRed'
_cm_GateRdBu1_data = (
        (0.0, _br), (0.5, _br), (0.5, _bb), (1.0, _bb)
        )
_cm_GateRdBu2_data = (
        (0.0, _br), (0.95, _br), (0.95, _bb), (1.0, _bb)
        )
_cm_GateBuRd2_data = (
        (0.0, _bb), (0.95, _bb), (0.95, _br), (1.0, _br)
        )
_cm_GateBuRd_data = (
        (0.0, 'Navy'), (0.5, _bb), (0.5, _br), (1.0, 'OrangeRed')
        )
_cm_GateBWR_data = (
        (0.0, '#0000FF'), (0.25, '#AFAFFF'), (0.25, '#144D8B'), (0.5, '#FFFFFF'), (0.75, '#B31B2C'), (0.75, '#FFAFAF'), (1.0, '#FF0000')
        )
_cm_GateBWR2_data = (
        (0.0, '#0000FF'), (0.1, '#AFAFFF'), (0.1, '#144D8B'), (0.5, '#FFFFFF'), (0.9, '#B31B2C'), (0.9, '#FFAFAF'), (1.0, '#FF0000')
        )
_cm_GateGBWRP_data = (
        (0.0, '#BBE28B'), (0.25, '#397B1C'), (0.25, '#0000FF'), (0.5, '#FFFFFF'), (0.75, '#FF0000'), (0.75, '#C1197A'), (1.0, '#FCDBEC')
        )
_cm_GateGBWRP2_data = (
        (0.0, '#BBE28B'), (0.1, '#397B1C'), (0.1, '#0000FF'), (0.5, '#FFFFFF'), (0.9, '#FF0000'), (0.9, '#C1197A'), (1.0, '#FCDBEC')
        )
_cm_pop_poster_data = list(pop_poster)
_cm_crayon_data = crayon[:14]
_cm_crayon_dark_data = crayon_dark[:10]

_cm_solid_data = (
        (0.0, '#A506A2'), (1.0/7.0, '#3A1FBE'), (2.0/7.0, '#045990'), (2.5/7.0, '#018950'),
        (3.0/7.0, '#24CC00'), (3.5/7.0, '#8DF200') , (4.0/7.0, '#FFFF00'),  (5.0/7.0, '#FFB500'), (6.0/7.0, '#FF6800'), 
        (6.5/7.0, '#FF0000'), (1.0, '#DF0000')
        )

_cm_night_data = [(0.0, '#27013a'), (0.1667, '#4763f7'), (0.3, '#82d6fa'), (0.45, '#92fee5'), (0.5, '#8ffe41'), (0.6667, '#efd743'), (0.85, '#cc2338'), (0.89, '#cd243e'), (0.94, '#e07ef0'), (0.985, '#f3d9f6'), (1.0, '#eedadc')]     
_cm_delight_data = ['#FFFFFF', '#874ac8', '#5419ad', '#2a36d6', '#0251fc', '#00e8ff', '#00f300', '#009000', '#5fc400', '#ffff00', '#ff7e00', '#ff0000', '#6d0000']            

_cm_halfgray1_data = ( (0.0, '1.0'), (1.0, '0.5') )
_cm_halfgray2_data = ( (0.0, '1.0'), (1.0, '0.3') )
_cm_halfgray3_data = ( (0.0, '1.0'), (0.35, '0.9'), (0.5, '0.8'), (0.7, '0.65'), (1.0, '0.3') )

todo_dict = {
        'cm_night':_cm_night_data,
        'cm_delight':_cm_delight_data,
        'cm_GateRdBu1':_cm_GateRdBu1_data,
        'cm_GateRdBu2':_cm_GateRdBu2_data,
        'cm_GateBuRd2':_cm_GateBuRd2_data,
        'cm_GateBuRd':_cm_GateBuRd_data,
        'cm_GateBWR':_cm_GateBWR_data,
        'cm_GateBWR2':_cm_GateBWR2_data,
        'cm_GateGBWRP':_cm_GateGBWRP_data,
        'cm_GateGBWRP2':_cm_GateGBWRP2_data,
        'cm_pop_poster':_cm_pop_poster_data,
        'cm_crayon':_cm_crayon_data,
        'cm_crayon_dark':_cm_crayon_dark_data,
        'cm_solid':_cm_solid_data,
        'cm_brewer_PuOr': brewer_PuOr,
        'cm_brewer_BrBG': brewer_BrBG,
        'cm_brewer_PRGn': brewer_PRGn,
        'cm_brewer_PiYG': brewer_PiYG,
        'cm_brewer_RdBu': brewer_RdBu,
        'cm_brewer_RdGy': brewer_RdGy,
        'cm_brewer_RdYlBu': brewer_RdYlBu,
        'cm_brewer_RdYlGn': brewer_RdYlGn,
        'cm_brewer_spectral': brewer_spectral,
        'cm_halfgray1': _cm_halfgray1_data,
        'cm_halfgray2': _cm_halfgray2_data,
        'cm_halfgray3': _cm_halfgray3_data,
        }

for _cm_name, _cm_data in todo_dict.items():
    if not _cm_name.endswith('_r') and _cm_name+'_r' not in todo_dict:
        if isinstance(_cm_data[0], tuple) and len(_cm_data[0]) == 2:
            # e.g. (0.0, (1.0, 0.0, 1.0))
            # or   (0.4, 'red')
            reversed_data = [(1.0-x, c) for x, c in reversed(_cm_data)]
        else:
            reversed_data = list(reversed(_cm_data))
    todo_dict[_cm_name + '_r'] = reversed_data

res_dict = dict()
for _cm_name, _cm_data in todo_dict.iteritems():
#    print _cm_name, _cm_data
    res_dict[_cm_name] = LSCM.from_list(_cm_name, _cm_data)

locals().update(res_dict)


__all__ = ['LSCM'] + [m for m in res_dict if m.startswith('cm_')]

if __name__ == '__main__':
    from matplotlib import pyplot as plt
    import numpy as np
    a=np.outer(np.arange(0,1,0.01),np.ones(10))
    plt.figure(figsize=(10,7))
    plt.subplots_adjust(top=0.75,bottom=0.05,left=0.01,right=0.99)
    maps=[m for m in todo_dict if not m.endswith('_r')]
    maps.sort()
    l=len(maps)+1
    for i, m in enumerate(maps):
        cmap=res_dict[m]
        plt.subplot(1,l,i+1)
        plt.axis("off")
        plt.imshow(a,aspect='auto',cmap=cmap,origin="lower")
        plt.title(m,rotation=90,fontsize=10, va='bottom')
    plt.savefig("metlib.color.cm.png",dpi=100,facecolor='w')
