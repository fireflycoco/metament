#!/usr/bin/env python
import numpy as np
from random import choice
from matplotlib import colors

def n_colors_towards_somecolor(basecolor, tocolor, N=5, max_ratio=0.7, use_basecolor_alpha=True):
    """Generate a list of N colors from basecolor to tocolor
    basecolor: starting color
    tocolor: to this color
    N: number of result colors
    max_ratio: 0-1, the degree of tranformation. if max_ratio == 0.0, all colors will be the same as the first color. if max_ratio == 1.0, the last color will be tocolor
    use_basecolor_alpha: use basecolor's alpha as all color's alpha
    """
    w_ratios = np.linspace(0.0, max_ratio, N)
    b_ratios = 1.0 - w_ratios
    converter = colors.ColorConverter()
    bcolor_rgba = np.array(converter.to_rgba(basecolor))
    tocolor_rgba = np.array(converter.to_rgba(tocolor))
    if use_basecolor_alpha:
        tocolor_rgba[-1] = bcolor_rgba[-1]
    res = []
    for i in range(N):
        res.append(tuple(bcolor_rgba*b_ratios[i] + tocolor_rgba*w_ratios[i]))
    return res

def n_colors_towards_white(basecolor, N=5, max_ratio=0.7):
    """Generate a list of N colors from basecolor to white
    basecolor: starting color
    N: number of result colors
    max_ratio: 0-1, the degree of tranformation. if max_ratio == 0.0, all colors will be the same as the first color. if max_ratio == 1.0, the last color will be tocolor
    """
    return n_colors_towards_somecolor(basecolor, 'white', N=N, max_ratio=max_ratio, use_basecolor_alpha=True)

class ColorBox(object):
    """ColorBox is a container for color list, which gives colors one by one or randomly when calling get_color() or __call__().
    """
    def __init__(self, colors):
        self.colors = colors
        self.i = 0

    def __call__(self, random=False):
        return self.get_color(random)

    def get_color(self, random=False):
        if random:
            return choice(self.colors)
        else:
            res = self.colors[self.i]
            self.i = (self.i + 1) % len(self.colors)
            return res

color_lists = ['reds', 'oranges', 'yellows', 'greens', 'blues', 'grays', 'greys', 'whites', 'purples', 'pinks', 'browns', 
        'seasons1', 'seasons2', 'seasons3', 'months1', 'months2', 'months3',
        'crayon', 'crayon_dark', 'crayon_dark_random',
        'pop_poster', 'pop_poster_random',
        'brewer_set1', 'brewer_set2', 'brewer_set3', 
        'brewer_pastel1', 'brewer_pastel2',
        'brewer_accent', 'brewer_dark2',
        'brewer_paired', 'brewer_paired2',
        'brewer_BrBG', 'brewer_PiYG', 'brewer_PRGn', 'brewer_PuOr',
        'brewer_RdBu', 'brewer_RdGy', 'brewer_RdYlBu', 'brewer_RdYlGn',
        'brewer_spectral'
        ]

reds = ['red', 'darkred', 'lightcoral', 'crimson', 'lightsalmon']
oranges = ['orange', 'orangered', 'sandybrown', 'goldenrod', 'wheat']
yellows = ['yellow', 'goldenrod', 'khaki', 'darkkhaki', 'moccasin']
greens = ['green', 'darkgreen', 'yellowgreen', 'olive', 'lightgreen', 'lime']
blues = ['blue', 'navy', 'royalblue', 'steelblue', 'skyblue', 'aqua']
grays = ['gray', 'dimgray', 'lightgrey', 'darkslategray', 'slategray', 'black']
greys = grays
whites = ['white', 'mistyrose', 'honeydew', 'antiquewhite', 'beige']
purples = ['purple', 'indigo', 'plum', 'blueviolet', 'mediumpurple']
pinks = ['hotpink', 'mediumvioletred', 'pink', 'magenta', 'lightcoral']
browns = ['sienna', 'saddlebrown', 'burlywood', 'brown', 'chocolate']

seasons1 = ['#006500', '#FF0F0F', '#DEB402', '#3E90FF']
seasons2 = ['darkgreen', 'red', 'goldenrod', 'blue']
seasons3 = ['lightgreen', 'salmon',  'khaki', 'deepskyblue']

months1 = []
months2 = []
for c in seasons2:
    months2.extend(n_colors_towards_white(c, 3, max_ratio=0.6))
months2 = months2[10:12] + months2[0:10]
months3 = []
for c in seasons3:
    months3.extend(n_colors_towards_somecolor(c, 'black', 3, max_ratio=0.5))
months3 = months3[10:12] + months3[0:10]

crayon = ['#EE204D', '#FFAACC', '#926EAE', '#1F75FE', '#1A4876', '#1CAC78', '#9FE2BF', '#BACF6C', '#FCE883', '#E7C697', '#FF7538', '#B4674D', '#EA7E5D', '#CB4154', '#95918C', 'black']
crayon_dark = ['#EE204D', '#926EAE', '#1F75FE', '#1A4876', '#1CAC78', '#BACF6C', '#FF7538', '#B4674D', '#EA7E5D', '#CB4154', '#95918C', '#202020']
crayon_dark_random = ['#1F75FE', '#CB4154', '#1CAC78', '#FF7538', '#926EAE', '#EA7E5D', '#202020', '#1A4876', '#EE204D', '#BACF6C', '#B4674D', '#95918C']

pop_poster = ['#861D2B', '#CA3226', '#E66B1E', '#FBD21A', '#BCC03E', '#69A646', '#318266', '#34A883', '#1F6597', '#273F61', '#8B4AA0', '#95879E', '#6E5248']
pop_poster_random = ['#1F6597', '#CA3226', '#69A646', '#E66B1E', '#8B4AA0', '#6E5248', '#34A883', '#273F61', '#861D2B', '#318266', '#95879E', '#FBD21A', '#BCC03E']

brewer_set1 = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999']        
brewer_set2 = [ '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3']        
brewer_set3 = ['#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F']
brewer_pastel1 = ['#FBB4AE', '#B3CDE3', '#CCEBC5', '#DECBE4', '#FED9A6', '#FFFFCC', '#E5D88D', '#FDDAEC', '#F2F2F2'] 
brewer_pastel2 = [ '#b3e2cd', '#fdcdac', '#cbd5e8', '#f4cae4', '#e6f5c9', '#fff2ae', '#f1e2cc', '#cccccc']        
brewer_accent = [ '#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f', '#bf5b17', '#666666']        
brewer_dark2 = [ '#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666']        
brewer_paired = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#BF5B17']
brewer_paired2 = ['#1F78B4', '#A6CEE3', '#33A02C', '#B2DF8A', '#E31A1C', '#FB9A99', '#FF7F00', '#FDBF6F', '#6A3D9A', '#CAB2D6', '#BF5B17', '#FFFF99']

brewer_BrBG = ['#543005', '#8c510a', '#bf812d', '#dfc27d', '#f6e8c3', '#f5f5f5', '#c7eae5', '#80cdc1', '#35978f', '#01665e', '#003c30']
brewer_PiYG = [ '#8e0152', '#c51b7d', '#de77ae', '#f1b6da', '#fde0ef', '#f7f7f7', '#e6f5d0', '#b8e186', '#7fbc41', '#4d9221', '#276419']
brewer_PRGn = [ '#40004b', '#762a83', '#9970ab', '#c2a5cf', '#e7d4e8', '#f7f7f7', '#d9f0d3', '#a6dba0', '#5aae61', '#1b7837', '#00441b']
brewer_PuOr = [ '#7f3b08', '#b35806', '#e08214', '#fdb863', '#fee0b6', '#f7f7f7', '#d8daeb', '#b2abd2', '#8073ac', '#542788', '#2d004b']
brewer_RdBu = [ '#67001f', '#b2182b', '#d6604d', '#f4a582', '#fddbc7', '#f7f7f7', '#d1e5f0', '#92c5de', '#4393c3', '#2166ac', '#053061']
brewer_RdGy = [ '#67001f', '#b2182b', '#d6604d', '#f4a582', '#fddbc7', '#ffffff', '#e0e0e0', '#bababa', '#878787', '#4d4d4d', '#1a1a1a']
brewer_RdYlBu = [ '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
brewer_RdYlGn = [ '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee08b', '#ffffbf', '#d9ef8b', '#a6d96a', '#66bd63', '#1a9850', '#006837']
brewer_spectral = [ '#9e0142', '#d53e4f', '#f46d43', '#fdae61', '#fee08b', '#ffffbf', '#e6f598', '#abdda4', '#66c2a5', '#3288bd', '#5e4fa2']
    

if __name__ == '__main__':
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from metlib.plot.axgrid import axgrid

    def plot_color_list(clist, figname=None):
        plt.figure()
        cmap = ListedColormap(clist, name='tmp_cmap')
        data = np.arange(len(clist)).reshape((1,len(clist)))
        plt.pcolor(data, cmap=cmap)
        if figname is None:
            plt.show()
        else:
            plt.savefig(figname)

    def plot_all_color_list(figname=None):
        fig = plt.figure(figsize=(8, 12))
        N = len(color_lists)
        axs = axgrid(fig, nrow=N, ncol=1,
                left=0.3, right=0.05, top=0.02, bottom=0.02,
                vspace=0.3/N).flatten()
        for ax in axs:
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
            plt.setp(ax.get_xticklines(), visible=False)
            plt.setp(ax.get_yticklines(), visible=False)
        for i, ax in enumerate(axs):
            cl = eval(color_lists[i])
            for ii, c in enumerate(cl):
                ax.axvspan(ii, ii+1, color=c)
            ax.set_ylabel(color_lists[i], fontsize='small', fontweight='semibold', 
                    rotation=0.0, ha='right', va='center')
            ax.set_xlim(0, len(cl))

        if figname is None:
            plt.show()
        else:
            plt.savefig(figname)
        

    if len(sys.argv) >= 3:
        figname = sys.argv[2]
    else:
        figname = None

    if len(sys.argv) >= 2:
        if sys.argv[1] == 'all':
            plot_all_color_list(figname)
        else:
            plot_color_list(sys.argv[1], figname)

