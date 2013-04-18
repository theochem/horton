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

from horton.cache import just_once
from horton.grid.cext import becke_helper_atom
from horton.log import log
from horton.part.base import DPart
from horton.periodic import periodic


__all__ = ['BeckeDPart']


class BeckeDPart(DPart):
    name = 'becke'
    options = ['local', 'k']

    '''Class for Becke partitioning'''
    def __init__(self, system, grid, local=True, k=3):
        self._k = k
        DPart.__init__(self, system, grid, local)

    def _init_log(self):
        DPart._init_log(self)
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Becke'),
                ('Switching function', 'k=%i' % self._k),
            ])
            log.cite('becke1988_multicenter', 'the use of Becke partitioning')

    @just_once
    def _init_partitioning(self):
        radii = np.array([periodic[n].cov_radius for n in self.system.numbers])
        for i in xrange(self.natom):
            grid = self.get_grid(i)
            at_weights, new = self.cache.load('at_weights', i, alloc=grid.size)
            if new:
                at_weights[:] = 1
                becke_helper_atom(grid.points, at_weights, radii, self.system.coordinates, i, self._k)

    def _get_k(self):
        '''The order of the Becke switching function.'''
        return self._k

    k = property(_get_k)
