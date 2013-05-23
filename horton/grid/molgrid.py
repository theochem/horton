# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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
'''Molecular integration grids'''



import numpy as np

from horton.grid.base import IntGrid
from horton.grid.atgrid import AtomicGrid, interpret_atspec
from horton.grid.cext import becke_helper_atom
from horton.log import log
from horton.periodic import periodic


__all__ = [
    'BeckeMolGrid'
]



class BeckeMolGrid(IntGrid):
    '''Molecular integration grid using Becke weights'''
    def __init__(self, system, atspecs='tv-13.1-3', k=3, random_rotate=True, keep_subgrids=False):
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
                By default, no (atomic) subgrids are stored. When set to True,
                atomic subgrids will be stored. This option is mainly of
                interest for AIM analysis.

           The argument atspecs may have two formats:

           * A single atspec tuple: ``(rgrid, nll)`` or ``(rgrid, nlls)``, where

             * ``rgrid`` is an instance of the RadialGrid
             * ``nll`` is the number of Lebedev-Laikov points on each sphere,
             * ``nlls`` is a list of numbers of Lebedev-Laikov grid points for
               the respective spheres of the atomic grid.

           * A list of atspec tuples as discussed in the foregoing bullet point,
             one for each atom. The length of this list must equal the number of
             atoms.
        '''
        # transform atspecs into usable format
        size, atspecs = get_mol_grid_size(system.numbers, atspecs, system.natom)

        # assign attributes
        self._system = system
        self._atspecs = atspecs
        self._keep_subgrids = keep_subgrids
        self._random_rotate = random_rotate
        self._k = k

        # allocate memory for the grid
        points = np.zeros((size, 3), float)
        weights = np.zeros(size, float)
        log.mem.announce(points.nbytes + weights.nbytes)

        # construct the atomic grids
        if keep_subgrids:
            atgrids = []
        else:
            atgrids = None
        offset = 0
        # More recent covalent radii are used than in the original work of Becke.
        cov_radii = np.array([periodic[n].cov_radius for n in system.numbers])
        for i in xrange(system.natom):
            rgrid, atnlls = atspecs[i]
            atsize = atnlls.sum()
            atgrid = AtomicGrid(system.numbers[i], system.coordinates[i],
                                (rgrid, atnlls), random_rotate,
                                points[offset:offset+atsize])
            weights[offset:offset+atsize] = atgrid.weights
            becke_helper_atom(points[offset:offset+atsize], weights[offset:offset+atsize], cov_radii, system.coordinates, i, self._k)
            if keep_subgrids:
                atgrids.append(atgrid)
            offset += atsize

        # finish
        IntGrid.__init__(self, points, weights, atgrids)

        # Some screen info
        self._log_init()

    def __del__(self):
        if log is not None and hasattr(self, 'weights'):
            log.mem.denounce(self.points.nbytes + self.weights.nbytes)

    def _get_system(self):
        '''The system object for which this grid is made.'''
        return self._system

    system = property(_get_system)

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

    def _log_init(self):
        if log.do_medium:
            log('Initialized: %s' % self)
            log.deflist([
                ('Size', self.size),
                ('System', self._system),
                ('Switching function', 'k=%i' % self._k),
            ])
            # Cite reference
            log.cite('becke1988_multicenter', 'the multicenter integration scheme used for the molecular integration grid')
            log.cite('cordero2008', 'the covalent radii used for the Becke-Lebedev molecular integration grid')



def get_mol_grid_size(numbers, atspecs, natom):
    '''Compute the size of the molecule grid and recreate atspecs with extra info.

       **Arguments:**

       numbers
            An array with element numbers. These are only used when atspecs
            is or contains string that refer to built-in grids.

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
        if not hasattr(atspecs[0], '__len__') or isinstance(atspecs, basestring):
            atspecs = [atspecs]*natom
        else:
            if not len(atspecs) == natom:
                raise TypeError('When a atspec is given per atom, the size of the list must match the number of atoms.')

    size = 0
    for i in xrange(len(atspecs)):
        rgrid, nlls = interpret_atspec(numbers[i], atspecs[i])
        atspecs[i] = rgrid, nlls
        size += nlls.sum()

    return size, atspecs
