#!/usr/bin/env python2.7

# layer.py

import os, sys, re
import numpy as np
from copy import copy, deepcopy
from .peak import *
from metlib.misc.datatype import Null, limited_int


class LayerCell(object):
    regex_str = r'Cell\[([^\[\]]*?)\]'
    regex = re.compile(regex_str)
    def __init__(self, ID, center, lower, upper, 
            sign = -1,
            life=10, fulllife=None, 
            healthy=True,
            peak=None,
            index=-1
            ):
        if fulllife is None:
            fulllife = life
        self.ID = ID
        self.center = center
        self.lower = lower
        self.upper = upper
        self.sign = sign
        self.life = limited_int(life, 0, fulllife)
        self.fulllife = fulllife
        self.peak = peak
        self.healthy = healthy
        self.index = index

    def __str__(self):
        signstr = '-' if self.sign == -1 else '+'
        healthystr = 'H' if self.healthy else 'U'
        return '[%s]{%s %s %s}%s%s%s' % (
            self.ID, self.lower, self.center, self.upper, 
            signstr, int(self.life), healthystr
            )

    def __repr__(self):
        signstr = '-' if self.sign == -1 else '+'
        healthystr = 'H' if self.healthy else 'U'
        return 'Cell[%d %d %d %d %d %s %d %d %s %s]' % (
                self.ID, self.index, 
                self.center, self.lower, self.upper,
                signstr, self.life, self.fulllife, healthystr,
                peak2str(self.peak) )


def _matchstr2cell(mstr):
    if len(mstr) == 0:
        cell = None
    else:
        tks = mstr.split(' ')
        try:
            sign = -1 if tks[5] == '-' else 1
            healthy = True if tks[8] == 'H' else False
            cell = LayerCell(ID=int(tks[0]), 
                    center=int(tks[2]), lower=int(tks[3]), upper=int(tks[4]),
                    sign=sign, life=int(tks[6]), fulllife=int(tks[7]),
                    healthy=healthy,
                    peak=parse_peak(tks[9]),
                    index=int(tks[1]))
        except Exception as e:
            cell = None
    return cell        

def parse_cell(s):
    m = re.search(LayerCell.regex, s)
    if m:
        return _matchstr2cell(m.group(1))
    else:
        return None

def parse_cells(s):
    res = []
    mstrs = re.findall(LayerCell.regex, s)
    for mstr in mstrs:
        res.append(_matchstr2cell(mstr))
    return res

def parse_cellfile(cell_fname):
    res = []
    with open(cell_fname) as cell_file:
        for line in cell_file:
            if line.lstrip().startswith('#'):
                continue
            cells = parse_cells(line)
            res.append(cells)
    return res

def cell2str(cell):
    if cell in (None, Null):
        return 'Cell[]'
    else:
        return repr(cell)

class LayerInfo(object):
    def __init__(self, beg_cell=None):
        if beg_cell is None:
            self.ID = None
            self.cells = []
        else:
            self.ID = beg_cell.ID
            self.cells = [beg_cell]
        self.attr = dict()
    
    # beg/end cells and index
    @property
    def beg_cell(self):
        return self.cells[0] if len(self) > 0 else None

    @property
    def end_cell(self):
        return self.cells[-1] if len(self) > 0 else None

    @property
    def beg_index(self):
        return self.cells[0].index if len(self) > 0 else None

    @property
    def end_index(self):
        return self.cells[-1].index if len(self) > 0 else None

    def append(self, cell):
        self.cells.append(cell)
    
    def __iter__(self):
        return self.cells.__iter__()

    def __getitem__(self, index):
        if isinstance(index, slice):
            oldstart=index.start
            oldstop = index.stop
            oldstep = index.step
            if oldstart is None:
                oldstart = 0
            if oldstep is None:
                oldstep = 1
            rel_index = slice(oldstart-self.beg_index, 
                    oldstop-self.beg_index,
                    oldstep)
            return self.cells[rel_index]
        rel_i = index - self.beg_index
        if 0 <= rel_i < len(self):
            return self.cells[rel_i]
        else:
            raise IndexError('Index %d out of layer range %d - %d' % \
                (index, self.beg_index, self.end_index))

    def show_layer(self):
        for cell in self.cells:
            print "$ %4d $" % cell.index, cell

    def __len__(self):
        return len(self.cells)

    def __str__(self):
        return "Layer %4d : {%4d %4d}" % (self.ID,  \
                self.beg_index, self.end_index)

    def __repr__(self):
        return self.__str__()

def cmp_layer(layer1, layer2):
    l1_info = dict([(i, c) for (i, c) in layer1])
    l2_info = dict([(i, c) for (i, c) in layer2])
    common_indexs = set(l1_info.keys()) & set(l2_info.keys())
    if len(common_indexs) == 0:
        return 0
    else:
        sample_index = common_indexs.pop()
        l1_cell = l1_info[sample_index]
        l2_cell = l2_info[sample_index]
        return cmp(l1_cell.center, l2_cell.center)

def get_layer_segment(layer, start_cell, max_number):
    start_index = start_cell.index
    last_index = start_index + max_number
    if last_index < layer.beg_index - 1:
        last_index = layer.beg_index -1
    step = int(np.sign(max_number))
    res = layer[start_index:last_index:step]
    if step == -1:
        res.reverse()
    return res

def intersect_width(lower1, upper1, lower2, upper2):
    if lower1 > upper2 or upper1 < lower2:
        res = 0
    else:
        u = min(upper1, upper2)
        l = max(lower1, lower2)
        res = u - l
    return res

class LayerMarker(object):
    regex_str = r'Marker\((.*?)\)'
    regex = re.compile(regex_str)

    def __init__(self, ID, layerID1, layerID2, x, y, marker):
        """marker: { } @"""
        self.ID = ID
        self.layerID1 = layerID1
        self.layerID2 = layerID2
        self.x = x
        self.y = y
        self.marker = marker

    def __cmp__(self, other):
        return cmp(self.x, other.x)

    def __repr__(self):
        return "Marker(%d,%d,%d,%d,%d,%s)" % (self.ID,
                self.layerID1, self.layerID2, 
                self.x, self.y, self.marker)
    
    @classmethod
    def parse(cls, s):
        try:
            match = re.search(cls.regex, s)
            tks = match.group(1).split(',')
            ID, ID1, ID2, x, y = [int(tk) for tk in tks[0:5]]
            marker = tks[5]
            return cls(ID, ID1, ID2, x, y, marker)
        except Exception as e:
            return None

    def remove(self):
        # for ArtistIDManager compatible
        pass

def save_layermarkers(fname, layermarkers):
    with open(fname, 'w') as outf:
        for m in sorted(layermarkers):
            if m is None:
                outf.write('Marker()')
            else:
                outf.write(repr(m))
            outf.write('\n')

def load_layermarkers(fname):
    res = set()
    with open(fname) as f:
        for line in f:
            l = line.strip()
            if len(line) == 0 or line.startswith('#'):
                continue
            res.add(LayerMarker.parse(l))
    return res

def markers2cells(guide_lc):
    markers = sorted(guide_lc.layermarkers)
    # add a leading '@' if not present.
    if len(markers) > 0 and markers[0].marker == '{' and markers[0].layerID1 != -1:
        ID = markers[0].layerID1
        layer = guide_lc.layers[ID]
        beg_cell = layer.beg_cell
        new_marker = LayerMarker(ID=0, layerID1=ID, layerID2=ID, 
                x=beg_cell.index, y=beg_cell.center, marker='@')
        markers.insert(0, new_marker)
    state = 'OFF'
    res_cells = []
    for i in range(len(markers)):
        m1 = markers[i]
        layerID = markers[i].layerID1
        if m1.marker in ('@', '}'):
            if (m1.marker == '@' and state == 'ON') or \
                    (m1.marker == '}' and m1.layerID1 == -1):
                state = 'OFF'
                continue
            beg_i = m1.x
            state = "ON"
            if i == len(markers) - 1:
                if layerID == -1:
                    end_i = len(guide_lc.indexs)
                else:
                    end_i = guide_lc.layers[layerID].end_index + 1
            else:
                m2 = markers[i+1]
                if m2.layerID1 != layerID:
                    end_i = min(guide_lc.layers[layerID].end_index + 1, m2.x + 1)
                    state = "OFF"
                elif m2.marker == '@':
                    end_i = m2.x + 1
                else:
                    end_i = m2.x
            sec_cells = guide_lc.layers[layerID][beg_i:end_i]
        elif markers[i].marker == '{':
            m1 = markers[i]
            m2 = markers[i+1]
            assert m2.marker == '}'
            # interp
            fill_xs = range(m1.x, m2.x)
            fill_centers = np.round(np.interp(fill_xs, [m1.x, m2.x], [m1.y, m2.y])).astype('i')
            sec_cells = [LayerCell(ID=-1, center=c, lower=c, upper=c,
                sign=-1, life=0, fulllife=10, healthy=False,
                peak=None, index=idx) for idx, c in zip(fill_xs, fill_centers) ]
            if m2.layerID1 == '-1':
                state = 'OFF'
            else:
                state = 'ON' 
        res_cells.extend(sec_cells)
    widths = [c.upper - c.lower for c in res_cells if c.ID != -1]
    mean_half_width = np.ma.masked_invalid(widths).mean() / 2.0
    if mean_half_width is np.ma.masked:
        mean_half_width = 0.0
    mean_half_width = np.round(mean_half_width).astype('i')
    for c in res_cells:
        if c.ID == -1:
            c.lower -= mean_half_width
            c.upper += mean_half_width
    return res_cells

