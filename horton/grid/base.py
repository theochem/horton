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
'''Base classes for 3D integration grids'''


import numpy as np

from horton.log import timer
from horton.grid.utils import parse_args_integrate
from horton.grid.cext import dot_multi, eval_spline_grid, \
    dot_multi_moments, eval_decomposition_grid
from horton.cext import Cell


__all__ = ['IntGrid']


class IntGrid(object):
    '''Base class for real-space integration grids in HORTON'''
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
        # assign begin and end attributes to the subgrids
        if subgrids is not None:
            offset = 0
            for sg in subgrids:
                sg.begin = offset
                offset += sg.size
                sg.end = offset

    def _get_size(self):
        '''The size of the grid.'''
        return self._weights.size

    size = property(_get_size)

    def _get_shape(self):
        '''The shape of the grid.'''
        return self._weights.shape

    shape = property(_get_shape)

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

    def zeros(self):
        return np.zeros(self.shape)

    def integrate(self, *args, **kwargs):
        '''Integrate the product of all arguments

           **Arguments:**

           data1, data2, ...
                All arguments must be arrays with the same size as the number
                of grid points. The arrays contain the functions, evaluated
                at the grid points, that must be multiplied and integrated.

           **Optional arguments:**

           center=None
                When given, multipole moments are computed with respect to
                this center instead of a plain integral.

           lmax=0
                The maximum angular momentum to consider when computing multipole
                moments

           mtype=1
                The type of multipole moments: 1=``cartesian``, 2=``pure``,
                3=``radial``, 4=``surface``.

           segments=None
                This argument can be used to divide the grid in segments. When
                given, it must be an array with the number of grid points in
                each consecutive segment. The integration is then carried out
                over each segment separately and an array of results is
                returned. The sum over all elements gives back the total
                integral.
        '''
        args, multipole_args, segments = parse_args_integrate(*args, **kwargs)
        args.append(self.weights)

        if multipole_args is None:
            # regular integration
            return dot_multi(*args, segments=segments)
        else:
            # computation of multipole expansion of the integrand
            center, lmax, mtype = multipole_args
            return dot_multi_moments(args, self.points, center, lmax, mtype, segments)


    @timer.with_section('Eval spher')
    def eval_spline(self, cubic_spline, center, output, cell=None):
        '''Evaluate a spherically symmetric function

           **Arguments:**

           cubic_spline
                A cubic spline with the radial dependence

           center
                The center of the spherically symmetric function

           output
                The output array

           **Optional arguments:**

           cell
                A unit cell when periodic boundary conditions are used.
        '''
        if cell is None:
            cell = Cell(None)
        eval_spline_grid(cubic_spline, center, output, self.points, cell)

    @timer.with_section('Eval decomp')
    def eval_decomposition(self, cubic_splines, center, output, cell=None):
        '''Evaluate a spherical decomposition

           **Arguments:**

           cubic_splines
                A list cubic splines, where each item is a radial function
                that is associated with a corresponding real spherical harmonic.

           center
                The center of the spherically symmetric function

           output
                The output array

           **Optional arguments:**

           cell
                A unit cell when periodic boundary conditions are used.
        '''
        if cell is None:
            cell = Cell(None)
        eval_decomposition_grid(cubic_splines, center, output, self.points, cell)
