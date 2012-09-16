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


from horton.grid.base import BaseGrid
from horton.grid.cext import lebedev_laikov_npoints
from horton.grid.sphere import LebedevLaikovSphereGrid


__all__ = [
    'AtomicGrid', 'get_atomic_grid_size',
]


class AtomicGrid(BaseGrid):
    def __init__(self, center, rtransform, int1d, nlls, random_rotate=True, points=None, keep_subgrids=0):
        '''
           **Arguments:**

           center
                The center of the radial grid

           rtransform
                An instance of a subclass of the BaseRTransform class.

           ind1d
                An instance of a subclass of the BaseIntegrator1D class.

           nlls
                The number Lebedev-Laikov grid points for each radial grid
                point. When this argument is not a list, all radial grid
                points get the same Lebedev-Laikov grid and the ``nsphere``
                argument must be given.

           **Optional arguments:**

           random_rotate
                When set to False, the random rotation of the grid points is
                disabled. Such random rotation improves the accuracy of the
                integration, but leads to small random changes in the results
                that are not reproducible.

           points
                Array to store the grid points

           keep_subgrids
                By default the (Lebedev-Laikov) subgrids are not stored
                separately. If set to 1, they are kept.
        '''
        size, nlls = get_atomic_grid_size(nlls, rtransform.npoint)
        self._center = center
        self._rtransform = rtransform
        self._int1d = int1d
        self._nlls = nlls
        self._random_rotate = random_rotate

        if points is None:
            points = np.zeros((size, 3), float)
        else:
            assert len(points) == size
        weights = np.zeros(size, float)

        if keep_subgrids > 0:
            llgrids = []
        else:
            llgrids = None
        offset = 0
        counter = 0
        nsphere = len(nlls)
        radii = rtransform.get_radii()
        rweights = int1d.get_weights(nsphere)
        rweights *= rtransform.get_volume_elements()
        for nll in nlls:
            llgrid = LebedevLaikovSphereGrid(center, radii[counter], nll, random_rotate, points[offset:offset+nll])
            if keep_subgrids > 0:
                llgrids.append(llgrid)
            weights[offset:offset+nll] = rweights[counter]*llgrid.weights
            offset += nll
            counter += 1

        BaseGrid.__init__(self, points, weights, llgrids)

    def _get_center(self):
        '''The center of the grid.'''
        return self._center

    center = property(_get_center)

    def _get_rtransform(self):
        '''The RTransform object of the grid.'''
        return self._rtransform

    rtransform = property(_get_rtransform)

    def _get_int1d(self):
        '''The 1D radial integrator object of the grid.'''
        return self._int1d

    int1d = property(_get_int1d)

    def _get_nlls(self):
        '''The number of Lebedev-Laikov grid points at each sphere.'''
        return self._nlls

    nlls = property(_get_nlls)

    def _get_nsphere(self):
        '''The number of spheres in the grid.'''
        return len(self._nlls)

    nsphere = property(_get_nsphere)

    def _get_random_rotate(self):
        '''The random rotation flag.'''
        return self._random_rotate

    random_rotate = property(_get_random_rotate)



def get_atomic_grid_size(nlls, nsphere):
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
        if len(nlls) != nsphere:
            raise ValueError('The size of the radial grid must match the number of elements in nlls')
    else:
        nlls = np.array([nlls]*nsphere, dtype=int)
    size = 0
    for nll in nlls:
        if nll not in lebedev_laikov_npoints:
            raise ValueError('A Lebedev-Laikov grid with %i points is not supported.' % nll)
        size += nll
    return size, nlls
