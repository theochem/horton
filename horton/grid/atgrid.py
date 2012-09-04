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

from horton.grid.cext import lebedev_laikov_npoints, lebedev_laikov_sphere


__all__ = [
    'AtomicGrid', 'get_atomic_grid_size', 'fill_atomic_grid',
    'get_random_rotation'
]


class AtomicGrid(object):
    def __init__(self, center, rtransform, nlls, nsphere=None, random_rotate=True):
        '''
           **Arguments:**

           center
                The center of the radial grid

           rtransform
                An instance of a subclass of the BaseRTransform class.

           nlls
                The number Lebedev-Laikov grid points for each radial grid
                point. When this argument is not a list, all radial grid
                points get the same Lebedev-Laikov grid and the nsphere
                argument must be given.

           **Optional arguments:**

           nsphere
                The number of radial grid points in the atomic grid, i.e. the
                number of spheres.

           random_rotate
                When set to False, the random rotation of the grid points is
                disabled. Such random rotation improves the accuracy of the
                integration, but leads to small random changes in the results
                that are not reproducible.
        '''
        size, nlls = get_atomic_grid_size(nlls, nsphere)
        self._center = center
        self._rtransform = rtransform
        self._nlls = nlls
        self._random_rotate = random_rotate

        self._points = np.zeros((size, 3), float)
        self._weights = np.zeros(size, float)
        fill_atomic_grid(self._points, self._weights, center, rtransform, nlls, random_rotate)

    def _get_center(self):
        '''The center of the grid.'''
        return self._center

    center = property(_get_center)

    def _get_rtransform(self):
        '''The RTransform object of the grid.'''
        return self._rtransform

    rtransform = property(_get_rtransform)

    def _get_nlls(self):
        '''The number of Lebedev-Laikov grid points at each sphere.'''
        return self._nlls

    nlls = property(_get_nlls)

    def _get_nsphere(self):
        '''The number of spheres in the grid.'''
        return len(self._nlls)

    nsphere = property(_get_nsphere)

    def _get_size(self):
        '''The size of the grid.'''
        return self._weights.shape[0]

    size = property(_get_size)

    def _get_random_rotate(self):
        '''The random rotation flag.'''
        return self._random_rotate

    random_rotate = property(_get_random_rotate)

    def _get_points(self):
        '''The grid points.'''
        return self._points

    points = property(_get_points)

    def _get_weights(self):
        '''The grid weights.'''
        return self._weights

    weights = property(_get_weights)



def get_atomic_grid_size(nlls, nsphere=None):
    '''Returns the size of the atomic grid.

       **Arguments:**

       nlls
            The number Lebedev-Laikov grid points for each radial grid
            point. When this argument is not a list, all radial grid
            points get the same Lebedev-Laikov grid and the nsphere
            argument must be given.

       **Optional arguments:**

       nsphere
            The number of radial grid points in the atomic grid, i.e. the
            number of spheres.

       **Returns:**

       size
            The total size of the atomic grid

       nlls
            A full list of Lebedev-Laikov grid sizes
    '''
    if hasattr(nlls, '__iter__'):
        nlls = np.array(nlls, dtype=int)
    else:
        if nsphere is None:
            raise ValueError('When nlls is not iterable, nsphere must be provided.')
        nlls = np.array([nlls]*nsphere, dtype=int)
    size = 0
    for nll in nlls:
        if nll not in lebedev_laikov_npoints:
            raise ValueError('A Lebedev-Laikov grid with %i points is not supported.')
        size += nll
    return size, nlls


def fill_atomic_grid(points, weights, center, rtransform, nlls, random_rotate):
    '''Fill the arrays points and weights with the atomic grid

       **Arguments:**

       points
            A numpy array with shape (size, 3), where size is obtained with
            the method ``get_atomic_grid_size``. This is an output array.

       weights
            A numpy array with shape (size,), where size is obtained with
            the method ``get_atomic_grid_size``. This is an output array.

       center
            The center of the radial grid

       rtransform
            An instance of a subclass of the BaseRTransform class.

       nlls
            The number Lebedev-Laikov grid points for each radial grid
            point. When this argument is not a list, all radial grid
            points get the same Lebedev-Laikov grid and the nsphere
            argument must be given.


       **Optional arguments:**

       random_rotate
            When set to False, the random rotation of the grid points is
            disabled. Such random rotation improves the accuracy of the
            integration, but leads to small random changes in the results
            that are not reproducible.
    '''
    offset = 0
    counter = 0
    nsphere = len(nlls)
    radii = rtransform.get_radii(nsphere)
    rweights = rtransform.get_int_weights(nsphere)
    for nll in nlls:
        lebedev_laikov_sphere(points[offset:offset+nll], weights[offset:offset+nll])
        if random_rotate:
            rotmat = get_random_rotation()
            points[offset:offset+nll] = np.dot(points[offset:offset+nll], rotmat)
        points[offset:offset+nll] *= radii[counter]
        weights[offset:offset+nll] *= rweights[counter]*4*np.pi*radii[counter]**2
        offset += nll
        counter += 1
    points += center


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
