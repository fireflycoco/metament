#!/usr/bin/env python2.7

# grid.py

import os, sys
import numpy as np
from matplotlib.patches import Polygon

class Grid2DBase(object):
    pass

class XYGrid2D(Grid2DBase):
    def __init__(self, X, Y):
        self._init(X, Y)

    def _init(self, X, Y):
        if np.ndim(X) == 1 and np.ndim(Y) == 1:
            self._X, self._Y = np.meshgrid(X, Y)
        elif np.ndim(X) == 1 and np.ndim(Y) == 2 and np.shape(X)[0] == np.shape(Y)[1]:
            self._X = np.tile(X, np.shape(Y)[0]).reshape(np.shape(Y))
            self._Y = Y
        elif np.ndim(X) == 2 and np.ndim(Y) == 1 and np.shape(X)[0] == np.shape(Y)[0]:
            self._X = X
            self._Y = np.repeat(Y, np.shape(X)[1]).reshape(np.shape(X))
        elif np.ndim(X) == 2 and np.ndim(Y) == 2 and np.shape(X) == np.shape(Y):
            self._X = X
            self._Y = Y
        else:
            msg = "X, Y shapes are not proper: X: %s, Y: %s" % (np.shape(X), np.shape(Y))
            raise ValueError(msg)

        corner_polygon = Polygon([
                                [self._X[0,0], self._Y[0,0]],
                                [self._X[0,-1], self._Y[0,-1]],
                                [self._X[-1,-1], self._Y[-1,-1]],
                                [self._X[-1,0], self._Y[-1,0]]
                                ])
        edge_polygon = Polygon(np.array(
            [np.concatenate((self._X[0,0:-1], self._X[0:-1,-1], self._X[-1,-1:0:-1], self._X[-1:0:-1,0])),
                np.concatenate((self._Y[0,0:-1], self._Y[0:-1,-1], self._Y[-1,-1:0:-1], self._Y[-1:0:-1,0]))]
                ).transpose())
        self._corner_path = corner_polygon.get_path()
        self._edge_path = edge_polygon.get_path()

    def corner_contains_points(self, points, radius=0.0):
        if np.shape(points) == (2,):
            return self._corner_path.contains_point(points, radius=radius)
        else:
            return self._corner_path.contains_points(points, radius=radius)

    def contains_points(self, points, radius=0.0):
        if np.shape(points) == (2,):
            return self._edge_path.contains_point(points, radius=radius)
        else:
            return self._edge_path.contains_points(points, radius=radius)

    @property
    def edge(self):
        return self._edge_path.vertices[:,0], self._edge_path.vertices[:,1]

    @property
    def X(self):
        return self._X

    @X.setter
    def X(self, new_X):
        self._init(new_X, self.Y)

    @property
    def Y(self):
        return self._Y

    @Y.setter
    def Y(self, new_Y):
        self._init(self.X, new_Y)

    def __str__(self):
        return "XYGrid2D:\nX:\n%s\nY:\n%s" % (str(self.X), str(self.Y))

    def __repr__(self):
        return self.__str__()

if __name__ == '__main__':
    pass
