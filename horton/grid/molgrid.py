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
'''Molecular integration grids'''



import numpy as np

from horton.grid.base import BaseGrid
from horton.grid.atgrid import AtomicGrid, get_atomic_grid_size
from horton.grid.cext import becke_helper_atom
from horton.periodic import periodic


__all__ = [
    'BeckeMolGrid'
]



class BeckeMolGrid(BaseGrid):
    '''Molecular integration grid using Becke weights'''
    def __init__(self, system, atspecs, k=3, random_rotate=True, keep_subgrids=0):
        '''
           **Arguments:**

           system
                The System object for which the molecular grid must be made.

           atspecs
                A specifications of the atomic grids. (See below.)

           **Optional arguments:**

           k
                The order of the switching function in Becke's weighting scheme.

           random_rotate
                Flag to control random rotation of spherical grids.

           keep_subgrids
                By default, not (atomic) subgrids are stored. When set to 1,
                atomic subgrids will be stored, but their (Lebedev-Laikov)
                subgrids are discarded. When set to 2, the Lebedev-Laikov
                subgrids are also stored internally. This option is mainly of
                interest for AIM analysis.

           The argument atspecs may have two formats:

           * A single atspec tuple: ``(rtransform, integrator1d, nll)``
             or ``(rtransform, integrator1d, nlls)``, where

             * ``rtransform`` is a transformation from a linear
               grid to some more convenient radial grid for integration and
             * ``integrator1d`` is an algorithm to integrate a function whose
               values are known on an equidistand grid with spacing 1.
             * ``nll`` is the number of Lebedev-Laikov points on each sphere,
             * ``nlls`` is a list of numbers of Lebedev-Laikov grid points for
               the respective spheres of the atomic grid.

           * A list of atspec tuples as discussed in the foregoing bullet point,
             one for each atom. The length of this list must equal the number of
             atoms.
        '''
        # transform atspecs into usable format
        size, atspecs = get_mol_grid_size(atspecs, system.natom)

        # assign attributes
        self._system = system
        self._atspecs = atspecs
        self._k = k
        self._random_rotate = random_rotate

        # allocate memory for the grid
        points = np.zeros((size, 3), float)
        weights = np.zeros(size, float)

        # construct the atomic grids
        if keep_subgrids > 0:
            atgrids = []
        else:
            atgrids = None
        offset = 0
        radii = np.array([periodic[n].cov_radius for n in system.numbers])
        for i in xrange(len(radii)):
            rtransform, int1d, atnlls, atsize = atspecs[i]
            atgrid = AtomicGrid(system.coordinates[i], rtransform, int1d,
                                atnlls, random_rotate,
                                points[offset:offset+atsize],
                                keep_subgrids=keep_subgrids-1)
            weights[offset:offset+atsize] = atgrid.weights
            becke_helper_atom(points[offset:offset+atsize], weights[offset:offset+atsize], radii, system.coordinates, i, self._k)
            if keep_subgrids > 0:
                atgrids.append(atgrid)
            offset += atsize

        # finish
        BaseGrid.__init__(self, points, weights, atgrids)

    def _get_system(self):
        '''The system object for which this grid is made.'''
        return self._system

    system = property(_get_system)

    def _get_coordinates(self):
        '''The coordinates of the grid.'''
        return self._coordinates

    coordinates = property(_get_coordinates)

    def _get_atspecs(self):
        '''The specifications of the atomic grids.'''
        return self._atspecs

    atspecs = property(_get_atspecs)

    def _get_k(self):
        '''The order of the Becke switching function.'''
        return self._k

    k = property(_get_k)

    def _get_random_rotate(self):
        '''The random rotation flag.'''
        return self._random_rotate

    random_rotate = property(_get_random_rotate)


def get_mol_grid_size(atspecs, natom):
    '''Compute the size of the molecule grid and recreate atspecs with extra info.

       **Arguments:**

       atspecs
            A list of specifications for the atomic grids.

       natom
            The total number of coordinates.
    '''
    if not hasattr(atspecs, '__len__'):
        raise TypeError('atspecs must be iterable')
    else:
        if len(atspecs) == 0:
            raise TypeError('atspecs must have at least one item.')
        if not hasattr(atspecs[0], '__len__'):
            atspecs = [atspecs]*natom
        else:
            if not len(atspecs) == natom:
                raise TypeError('When a atspec is given per atom, the size of the list must match the number of atoms.')

    size = 0
    for i in xrange(len(atspecs)):
        atspec = atspecs[i]
        if len(atspec) != 3:
            raise TypeError('An atomic grid spec must contain three elements.')
        rtransform, int1d, nlls = atspec
        if not hasattr(nlls, '__len__'):
            nlls = [nlls]*rtransform.npoint
        atsize = get_atomic_grid_size(nlls, rtransform.npoint)[0]
        atspecs[i] = rtransform, int1d, nlls, atsize
        size += atsize

    return size, atspecs
