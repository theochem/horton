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
from horton.dpart.base import BaseDPart
from horton.log import log


__all__ = ['HirshfeldDPart']


class HirshfeldDPart(BaseDPart):
    '''Base class for Hirshfeld partitioning'''
    def __init__(self, molgrid, proatomdb, local=True):
        self._proatomdb = proatomdb
        BaseDPart.__init__(self, molgrid, local)

    def _init_log(self):
        BaseDPart._init_log(self)
        if log.do_medium:
            with log.section('DPART'):
                log('  Scheme:            &Hirshfeld')
                log('  Proatomic DB:      &%s' % self._proatomdb)

    @just_once
    def _init_proatom_fns(self):
        for i in xrange(self.system.natom):
            proatom_fn = self.cache.load('proatom_fn', i, default=None)
            if proatom_fn is None:
                n = self.system.numbers[i]
                proatom_fn = self._proatomdb.get_hirshfeld_proatom_fn(n)
                self.cache.dump('proatom_fn', i, proatom_fn)

    @just_once
    def _init_at_weights(self):
        self._init_proatom_fns()
        # TODO: optimize for global grid
        for i0, grid in self.iter_grids():
            at_weights, new = self.cache.load('at_weights', i0, alloc=grid.size)
            if new:
                work = np.zeros(grid.size, float)
                pro_mol = np.zeros(grid.size, float)
                d = np.zeros(grid.size, float)
                for i1 in xrange(self.system.natom):
                    grid.distances(self.system.coordinates[i1], d)
                    proatom_fn = self.cache.load('proatom_fn', i1)
                    proatom_fn(d, work)
                    if i1 == i0:
                        at_weights[:] = work
                    pro_mol[:] += work
                # The following seems worse than it is. It does nothing to the
                # relevant numbers. It just avoids troubles in the division.
                pro_mol[:] += 1e-100
                at_weights[:] /= pro_mol
