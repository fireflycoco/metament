import numpy as np
from metlib.misc.datatype import isinteger

__all__ = ['ArtistIDManager']

class ArtistIDManager(object):
    """ArtistIDManager is a dict-like container which manages Artists in groups (ID).
Examples
--------
>>> m = ArtistIDManager()
>>> group_ID_1 = 1
>>> m[group_ID_1] = some_line      # add a line artist to group 1
>>> m[group_ID_1] = some_polygon   # add a polygon artist to the same group
>>> print m[group_ID_1]            # look up artists with ID
[some_line, some_polygon]
>>> print m[some_line]             # look up ID with artist
1
>>> group_ID_1 in m
True
>>> some_polygon in m
True
>>> del m[group_ID_1]              # del that group (and remove the artists from the fig).
    """
    
    def __init__(self):
        self.artist2ID = dict()
        self.ID2artist = dict()
        self.nowID = 1

    def isID(self, key):
        return isinstance(key, (int, long, np.integer, str, unicode, np.string_))

    def get_new_ID(self):
        res = self.nowID
        self.nowID += 1
        return res

    def __getitem__(self, key):
        if self.isID(key):
            return self.ID2artist[key]
        else:
            return self.artist2ID[key]

    def __setitem__(self, key, value):
        if self.isID(key):
            if key in self.ID2artist:
                self.ID2artist[key].append(value)
            else:
                self.ID2artist[key] = [value]
            self.artist2ID[value] = key
            if isinteger(key) and key >= self.nowID:
                self.nowID = key + 1
        else:
            self.__setitem__(value, key)

    def __delitem__(self, key):
        if not self.isID(key):
            if key not in self.artist2ID:
                return
            key = self.artist2ID[key]
        if key not in self.ID2artist:
            return
        for artist in self.ID2artist[key]:
            artist.remove()
            del self.artist2ID[artist]
        del self.ID2artist[key]

    def __contains__(self, key):
        return key in self.artist2ID or key in self.ID2artist

    def clear(self):
        for artist in self.artist2ID:
            artist.remove()
        self.artist2ID.clear()
        self.ID2artist.clear()
        self.nowID = 1

