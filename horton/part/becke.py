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
'''Becke partitioning'''


import numpy as np

from horton.grid.cext import becke_helper_atom
from horton.log import log, timer
from horton.part.base import WPart
from horton.periodic import periodic
from horton.units import angstrom


__all__ = ['BeckeWPart']


class BeckeWPart(WPart):
    '''Becke partitioning with Becke-Lebedev grids'''

    name = 'b'
    options = ['lmax', 'k']
    linear = True

    def __init__(self, coordinates, numbers, pseudo_numbers, grid, moldens,
                 spindens=None, local=True, lmax=3, k=3):
        '''
           **Optional arguments:** (that are not defined in ``WPart``)

           k
                The order of the polynomials used in the Becke partitioning.
        '''
        self._k = k
        WPart.__init__(self, coordinates, numbers, pseudo_numbers, grid,
                       moldens, spindens, local, lmax)

    def _init_log_scheme(self):
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Becke'),
                ('Switching function', 'k=%i' % self._k),
            ])
            log.cite('becke1988_multicenter', 'the use of Becke partitioning')
            log.cite('slater1964', 'the Brag-Slater radii used in the Becke partitioning')

    @timer.with_section('Becke part')
    def update_at_weights(self):
        if log.do_medium:
            log('Computing Becke weights.')

        # The list of radii is constructed to be as close as possible to
        # the original values used by Becke.
        radii = []
        for number in self.numbers:
            if number == 1:
                radius = 0.35*angstrom # exception defined in Becke's paper
            else:
                radius = periodic[number].becke_radius
                if radius is None: # for cases not covered by Brag-Slater
                    radius = periodic[number].cov_radius
            radii.append(radius)
        radii = np.array(radii)

        # Actual work
        pb = log.progress(self.natom)
        for index in xrange(self.natom):
            grid = self.get_grid(index)
            at_weights = self.cache.load('at_weights', index, alloc=grid.shape)[0]
            at_weights[:] = 1
            becke_helper_atom(grid.points, at_weights, radii, self.coordinates, index, self._k)
            pb()

    def _get_k(self):
        '''The order of the Becke switching function.'''
        return self._k

    k = property(_get_k)
