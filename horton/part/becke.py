# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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
"""Becke partitioning"""


from __future__ import print_function

import numpy as np

from .base import WPart
from .utils import angstrom, radius_becke, radius_covalent
from horton.grid import becke_helper_atom


__all__ = ['BeckeWPart']


class BeckeWPart(WPart):
    """Becke partitioning with Becke-Lebedev grids"""

    name = 'b'
    options = ['lmax', 'k']
    linear = True

    def __init__(self, coordinates, numbers, pseudo_numbers, grid, moldens,
                 spindens=None, local=True, lmax=3, k=3):
        """
           **Optional arguments:** (that are not defined in ``WPart``)

           k
                The order of the polynomials used in the Becke partitioning.
        """
        self._k = k
        WPart.__init__(self, coordinates, numbers, pseudo_numbers, grid,
                       moldens, spindens, local, lmax)

    def _init_log_scheme(self):
        print('5: Initialized: %s' % self)
        print([
            ('5: Scheme', 'Becke'),
            ('5: Switching function', 'k=%i' % self._k),
        ])
        self.biblio.append(['becke1988_multicenter', 'the use of Becke partitioning'])
        self.biblio.append(['slater1964', 'the Brag-Slater radii used in the Becke partitioning'])

    def update_at_weights(self):
        print('5:Computing Becke weights.')

        # The list of radii is constructed to be as close as possible to
        # the original values used by Becke.
        radii = []
        for number in self.numbers:
            if number == 1:
                # exception defined in Becke's paper
                radius = 0.35 * angstrom
            else:
                radius = radius_becke[number]
                if radius is None:
                    # for cases not covered by Brag-Slater
                    radius = radius_covalent[number]
            radii.append(radius)
        radii = np.array(radii)

        # Actual work
        for index in range(self.natom):
            grid = self.get_grid(index)
            at_weights = self.cache.load('at_weights', index, alloc=grid.shape)[0]
            at_weights[:] = 1
            becke_helper_atom(grid.points, at_weights, radii, self.coordinates, index, self._k)

    def _get_k(self):
        """The order of the Becke switching function."""
        return self._k

    k = property(_get_k)
