# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
'''Grids suitable for visualization'''


import numpy as np

from horton.grid.base import IntGrid


__all__ = ['LineGrid', 'RectangleGrid']


class LineGrid(IntGrid):
    '''Grid points evenly distributed over a 3D line segment'''
    def __init__(self, p1, p2, size, extend=0):
        '''
           **Arguments:**

           p1, p2
                Two points defining the line segment

           size
                The number of grid points

           **Optional arguments:**

           extend
                Used to control how far the grid points extrapolate from the
                line segment.
        '''
        self._p1 = p1
        self._p2 = p2
        norm = np.linalg.norm(p2 - p1)
        self._x = norm*(np.arange(size, dtype=float)/(size-1)-0.5)
        self._x *= 1 + 2*extend
        points = np.outer(self._x, (p2 - p1)/norm) + 0.5*(p2 + p1)
        assert points.shape == (size, 3)
        weight = norm/(size-1)
        weights = np.empty(size)
        weights.fill(weight)
        IntGrid.__init__(self, points, weights)

    def _get_p1(self):
        '''The first point of the line segment'''
        return self._p1.view()

    p1 = property(_get_p1)

    def _get_p2(self):
        '''The second point of the line segment'''
        return self._p2.view()

    p2 = property(_get_p2)

    def _get_x(self):
        '''The 1D axis for the grid points, useful for plotting

           :math:`p_1` and :math:`p_2` correspond to :math:`x_1=-|p_2-p_1|/2`
           and :math:`x_2=+|p_2-p_1|/2`, respectively.
        '''
        return self._x.view()

    x = property(_get_x)


class RectangleGrid(IntGrid):
    '''A 2D rectungular grid in 3D space'''
    def __init__(self, origin, axis0, axis1, l0, h0, l1, h1):
        '''
           **Arguments:**

           origin
                The origin for the definition of the rectangle.

           axis0, axis1
                The basis vectors that define the plane of the rectangle and
                the 2D space in this plane.

           l0, h0
                The lowest and highest position along axis0. (These must be
                integers.) Hence, along the first axis, there are (h0-l0+1)
                grid points.

           l1, h1
                The lowest and highest position along axis1.
        '''
        if h0 <= l0:
            raise ValueError('l0 should be lower than h0.')
        if h1 <= l1:
            raise ValueError('l1 should be lower than h1.')

        self._origin = origin
        self._axis0 = axis0
        self._axis1 = axis1
        self._l0 = l0
        self._h0 = h0
        self._l1 = l1
        self._h1 = h1

        size = (h0-l0+1)*(h1-l1+1)
        points = np.zeros((size, 3), float)
        weights = np.empty(size, float)
        weights[:] = np.sqrt(
            (np.linalg.norm(axis0)*np.linalg.norm(axis1))**2 -
            np.dot(axis0, axis1)**2
        )
        counter = 0
        for i0 in xrange(l0, h0+1):
            for i1 in xrange(l1, h1+1):
                points[counter] = origin + axis0*i0 + axis1*i1
                counter += 1
        IntGrid.__init__(self, points, weights)

    def _get_origin(self):
        '''The origin for the definition of the rectangle'''
        return self._origin.view()

    origin = property(_get_origin)

    def _get_axis0(self):
        '''The first basis vector'''
        return self._axis0.view()

    axis0 = property(_get_axis0)

    def _get_axis1(self):
        '''The second basis vector'''
        return self._axis1.view()

    axis1 = property(_get_axis1)

    def _get_l0(self):
        '''The lowest index along axis0'''
        return self._l0

    l0 = property(_get_l0)

    def _get_h0(self):
        '''The highest index along axis0'''
        return self._h0

    h0 = property(_get_h0)

    def _get_l1(self):
        '''The lowest index along axis1'''
        return self._l1

    l1 = property(_get_l1)

    def _get_h1(self):
        '''The highest index along axis1'''
        return self._h1

    h1 = property(_get_h1)

    def prepare_contour(self, data):
        '''Returns arguments suitable for matplotlib.pyplot.contour

           **Arguments:**

           data
                A vector with the correct number of elements (self.size)

           **Returns:** suitable arguments for  matplotlib.pyplot.contour

           x, y
                A numpy vector with the x and y axes (measured in units defined
                by the lengths of axis0 and axis1, respectively.

           z
                the same as data, transformed in a suitable 2D numpy array
        '''
        x = np.arange(self.l0, self.h0+1)*np.linalg.norm(self.axis0)
        y = np.arange(self.l1, self.h1+1)*np.linalg.norm(self.axis1)
        z = data.reshape((self.h0 - self.l0 + 1, self.h1 - self.l1 + 1))
        return x, y, z
