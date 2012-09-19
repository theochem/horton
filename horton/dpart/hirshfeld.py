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
from horton.grid.int1d import TrapezoidIntegrator1D
from horton.grid.cext import dot_multi
from horton.log import log


__all__ = ['HirshfeldDPart', 'HirshfeldIDPart']


class HirshfeldDPart(DPart):
    '''Base class for Hirshfeld partitioning'''
    def __init__(self, molgrid, proatomdb, local=True):
        self._proatomdb = proatomdb
        self._pro_mol_valid = False
        DPart.__init__(self, molgrid, local)

    def _init_log(self):
        DPart._init_log(self)
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld'),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('hirshfeld1977', 'for the use of Hirshfeld partitioning')

    def _at_weights_helper(self, i0, grid, at_weights):
        # Load work arrays from cache for efficiency
        (work, pro_mol, d), new = self.cache.load('at_weights_work', grid.size, alloc=(3, grid.size))
        if not new:
            pro_mol[:] = 0.0
        if self.local or not self._pro_mol_valid:
            # In case of local grids, the pro-molecule must be recomputed for
            # every grid. In case of a global grid, this only happens the
            # first time after the pro-atoms were updated.
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
        else:
            # In case of a global grid, and when the pro-molecule is up to date,
            # only the pro-atom needs to be recomputed.
            grid.distances(self.system.coordinates[i0], d)
            proatom_fn = self.cache.load('proatom_fn', i0)
            proatom_fn(d, work)
            at_weights[:] = work
        # Finally compute the ratio
        at_weights[:] /= pro_mol

    def _at_weights_cleanup(self):
        # Get rid of cached work arrays
        for i, grid in self.iter_grids():
            self.cache.discard('at_weights_work', id(grid.size))

    @just_once
    def _init_at_weights(self):
        for i in xrange(self.system.natom):
            proatom_fn = self.cache.load('proatom_fn', i, default=None)
            if proatom_fn is None:
                n = self.system.numbers[i]
                proatom_fn = self._proatomdb.get_hirshfeld_proatom_fn(n)
                self.cache.dump('proatom_fn', i, proatom_fn)

        for i0, grid in self.iter_grids():
            at_weights, new = self.cache.load('at_weights', i0, alloc=grid.size)
            if new:
                self._at_weights_helper(i0, grid, at_weights)

        self._at_weights_cleanup()


class HirshfeldIDPart(HirshfeldDPart):
    '''Iterative Hirshfeld partitioning'''
    def __init__(self, molgrid, proatomdb, local=True, threshold=1e-4, maxiter=500):
        self._threshold = threshold
        self._maxiter = maxiter
        HirshfeldDPart.__init__(self, molgrid, proatomdb, local)

    def _init_log(self):
        DPart._init_log(self)
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld-I'),
                ('Convergence threshold', '%.1e' % self._threshold),
                ('Maximum iterations', self._maxiter),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('bultinck2007', 'for the use of Hirshfeld-I partitioning')

    def _proatom_change(self, old, new):
        # Note that the dumb TrapezoidIntegrator is fine enough for
        # the convergence test.
        return (4*np.pi)*dot_multi(
            self._proatomdb._rtransform.get_radii()**2, # TODO: get routines are slow
            self._proatomdb._rtransform.get_volume_elements(),
            TrapezoidIntegrator1D().get_weights(self._proatomdb._rtransform.npoint),
            (old.copy_y() - new.copy_y())**2 # TODO: copy is slow
        )

    @just_once
    def _init_at_weights(self):
        # Perform one general check in the beginning to keep things simple.
        if all(self.cache.has('at_weights', i) for i in xrange(self.system.natom)):
            return

        self.do_mol_dens()

        # Iterative Hirshfeld loop
        populations = self.system.numbers.astype(float)
        counter = 0
        if log.do_medium:
            log('Iterative Hirshfeld partitioning loop')
            log.hline()
            log('Counter      Change   Pop.Error')
            log.hline()
        while True:
            # Update current pro-atoms
            change = 0.0
            first_iter = True
            for i in xrange(self.system.natom):
                old_proatom_fn = self.cache.load('proatom_fn', i, default=None)
                n = self.system.numbers[i]
                p = populations[i]
                proatom_fn = self._proatomdb.get_hirshfeld_i_proatom_fn(n, p)
                self.cache.dump('proatom_fn', i, proatom_fn)
                if old_proatom_fn is not None:
                    first_iter = False
                    change += self._proatom_change(old_proatom_fn, proatom_fn)

            # Enforce (single) update of pro-molecule in case of a global grid
            if not self.local:
                self._pro_mol_valid = False
            # Compute populations
            for i, grid in self.iter_grids():
                # Compute weight
                at_weights, new = self.cache.load('at_weights', i, alloc=grid.size)
                self._at_weights_helper(i, grid, at_weights)

                # Compute population
                dens = self.cache.load('mol_dens', i)
                populations[i] = grid.integrate(at_weights, dens)

            change = np.sqrt(change/self.system.natom)
            if log.do_medium:
                pop_error = populations.sum() - self.system.wfn.nel
                log('%7i  %10.5e  %10.5e' % (counter, change, pop_error))

            if not first_iter and change < self._threshold:
                self._converged = True
                break

            if counter > self._maxiter:
                break

            counter += 1

        if log.do_medium:
            log.hline()
            log('Converged: %s' % self._converged)
            log.blank()

        self._at_weights_cleanup()
        self.cache.dump('populations', populations)
