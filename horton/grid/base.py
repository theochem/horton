# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


from horton.grid.cext import dot_multi, grid_distances, eval_spline_grid
from horton.cext import Cell


__all__ = ['IntGrid']


class IntGrid(object):
    '''Base class for real-space integration grids in Horton'''
    def __init__(self, points, weights, subgrids=None):
        '''
           **Arguments:**

           points
                A numpy array with shape (npoint,3) with the Cartesian
                coordinates of the grids points.

           weights
                The integration weights of the grid points

           **Optional arguments:**

           subgrids
                Can be given when this grid is composed of several other
                grids. The points data is shared, but the weights of the
                subgrids may be different.
        '''
        if subgrids is not None and len(subgrids) == 0:
            raise TypeError('When subgrids are given, it may not be an empty list.')
        self._points = points
        self._weights = weights
        self._subgrids = subgrids

    def _get_size(self):
        '''The size of the grid.'''
        return self._weights.shape[0]

    size = property(_get_size)

    def _get_points(self):
        '''The grid points.'''
        return self._points

    points = property(_get_points)

    def _get_weights(self):
        '''The grid weights.'''
        return self._weights

    weights = property(_get_weights)

    def _get_subgrids(self):
        '''A list of grid objects used to construct this grid.'''
        return self._subgrids

    subgrids = property(_get_subgrids)

    def integrate(self, *args):
        '''Integrate the product of all arguments

           **Arguments:**

           data1, data2, ...
                All arguments must be arrays with the same size as the number
                of grid points. The arrays contain the functions, evaluated
                at the grid points, that must be multiplied and integrated.
        '''
        args = [arg for arg in args if arg is not None]
        return dot_multi(self.weights, *args)

    def distances(self, center, d):
        '''Compute distances between all grid points and a center, store result in d.'''
        # TODO: this should no longer be used. may still be useful for testing
        grid_distances(self.points, center, d)

    def eval_spline(self, cubic_spline, center, output, cell=None):
        if cell is None:
            cell = Cell(None)
        eval_spline_grid(cubic_spline, center, output, self.points, cell)
