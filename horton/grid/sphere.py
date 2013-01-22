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


import numpy as np

from horton.grid.base import IntGrid
from horton.grid.cext import lebedev_laikov_npoints, lebedev_laikov_sphere


__all__ = [
    'LebedevLaikovSphereGrid', 'get_random_rotation'
]


class LebedevLaikovSphereGrid(IntGrid):
    '''A spherical Lebedev-Laikov grid'''
    def __init__(self, center, radius, nll, random_rotate=True, points=None):
        '''
           **Arguments:**

           center
                The center of the sphere

           radius
                The radius of the sphere

           nll
                The number of Lebedev-Laikov grid points

           **Optional arguments:**

           random_rotate
                Flag to control random rotation of spherical grids.

           points
                Array to store the grid points
        '''
        self._center = center
        self._radius = radius
        self._nll = nll
        self._random_rotate = random_rotate

        if nll not in lebedev_laikov_npoints:
            raise ValueError('Unsupported number of Lebedev-Laikov grid points: %i' % nll)

        if points is None:
            points = np.zeros((nll, 3), float)
        else:
            assert len(points) == nll
        weights = np.zeros(nll, float)

        lebedev_laikov_sphere(points, weights)
        points *= radius
        if random_rotate:
            rotmat = get_random_rotation()
            points[:] = np.dot(points, rotmat)
        points[:] += center
        weights *= 4*np.pi*radius**2

        IntGrid.__init__(self, points, weights, None)

    def _get_center(self):
        '''The center of the sphere.'''
        return self._center

    center = property(_get_center)

    def _get_radius(self):
        '''The radius of the sphere.'''
        return self._radius

    radius = property(_get_radius)

    def _get_nll(self):
        '''The number of Lebedev-Laikov grid points.'''
        return self._nll

    nll = property(_get_nll)

    def _get_random_rotate(self):
        '''The random rotation flag.'''
        return self._random_rotate

    random_rotate = property(_get_random_rotate)



def get_random_rotation():
    '''Return a random rotation matrix'''
    # Get a random unit vector for the axis
    while True:
        axis = np.random.uniform(-1, 1, 3)
        norm = np.linalg.norm(axis)
        if norm < 1.0 and norm > 0.1:
            axis /= norm
            break
    x, y, z = axis

    # Get a random rotation angle
    angle = np.random.uniform(0, 2*np.pi)
    c = np.cos(angle)
    s = np.sin(angle)

    # Rodrigues' rotation formula
    return np.array([
        [x*x*(1-c)+c  , x*y*(1-c)-z*s, x*z*(1-c)+y*s],
        [x*y*(1-c)+z*s, y*y*(1-c)+c  , y*z*(1-c)-x*s],
        [x*z*(1-c)-y*s, y*z*(1-c)+x*s, z*z*(1-c)+c  ],
    ])
