# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import numpy as np
cimport numpy as np
np.import_array()

cimport lebedev_laikov

__all__ = [
    # lebedev_laikov
    'lebedev_laikov_npoint', 'lebedev_laikov_sphere', 'lebedev_laikov_npoints'
]


lebedev_laikov_npoints = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810]

def lebedev_laikov_npoint(int lvalue):
    '''lebedev_laikov_npoint(lvalue)

       Return the number of Lebedev-Laikov grid points for a given angular
       momentum.
    '''
    return lebedev_laikov.lebedev_laikov_npoint(lvalue)

def lebedev_laikov_sphere(np.ndarray[double, ndim=2] points,
                          np.ndarray[double, ndim=1] weights):
    '''lebedev_laikov_sphere(grid)

       Fill the grid with a Lebedev Laikov grid points of a given size.

       **Arguments:**

       points
            The output array for the grid points, shape (npoint,3).

       weights
            The output array for the grid weights, shape (npoint,).
    '''
    assert points.flags['C_CONTIGUOUS']
    assert points.shape[1] == 3
    assert points.shape[0] == weights.shape[0]
    lebedev_laikov.lebedev_laikov_sphere(points.shape[0], <double*>points.data,
                                         <double*>weights.data)
