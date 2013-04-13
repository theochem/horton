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
from horton.dpart.base import DPart
from horton.log import log


__all__ = ['HirshfeldDPart']


# TODO: proofread and add tests for pseudo densities


class HirshfeldDPart(DPart):
    name = 'h'
    options = ['local']

    '''Base class for Hirshfeld partitioning'''
    def __init__(self, molgrid, proatomdb, local=True):
        self._proatomdb = proatomdb
        DPart.__init__(self, molgrid, local)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def _init_log(self):
        DPart._init_log(self)
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld'),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('hirshfeld1977', 'the use of Hirshfeld partitioning')

    def _compute_at_weights(self, i0, grid, at_weights):
        promol, new = self.cache.load('promol', grid.size, alloc=grid.size)
        if new or self.local:
            # In case of local grids, the pro-molecule must always be recomputed.
            promol[:] = 0.0 # needed if not new and local.
            for i1 in xrange(self.system.natom):
                proatom_fn = self.cache.load('proatom_fn', i1)
                if i1 == i0:
                    at_weights[:] = 0.0
                    grid.eval_spline(proatom_fn, self.system.coordinates[i1], at_weights)
                    promol += at_weights
                else:
                    grid.eval_spline(proatom_fn, self.system.coordinates[i1], promol)
            # The following seems worse than it is. It does nothing to the
            # relevant numbers. It just avoids troubles in the division.
            promol[:] += 1e-100
        else:
            # In case of a global grid and when the pro-molecule is up to date,
            # only the pro-atom needs to be recomputed.
            proatom_fn = self.cache.load('proatom_fn', i0)
            at_weights[:] = 0.0
            grid.eval_spline(proatom_fn, self.system.coordinates[i0], at_weights)
        # Finally compute the ratio
        at_weights[:] /= promol

    def _at_weights_cleanup(self):
        # Get rid of cached work arrays
        for i, grid in self.iter_grids():
            self.cache.discard('at_weights_work', id(grid.size))

    @just_once
    def _init_at_weights(self):
        # TODO. Splines should not be stored in the cache. Use same structure as
        # in CPart: get_proatom_spline and _*_propars methods.
        # Create the splines and store them in the cache
        for i in xrange(self.system.natom):
            proatom_fn = self.cache.load('proatom_fn', i, default=None)
            if proatom_fn is None:
                n = self.system.numbers[i]
                proatom_fn = self._proatomdb.get_spline(n)
                self.cache.dump('proatom_fn', i, proatom_fn)

        # Compute the weights
        for i0, grid in self.iter_grids():
            at_weights, new = self.cache.load('at_weights', i0, alloc=grid.size)
            if new:
                self._compute_at_weights(i0, grid, at_weights)

        self._at_weights_cleanup()
