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
from horton.dpart.hirshfeld import HirshfeldDPart
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import dot_multi
from horton.log import log


__all__ = ['HirshfeldIDPart']


# TODO: proofread and add tests for pseudo densities


class HirshfeldIDPart(HirshfeldDPart):
    name = 'hi'
    options = ['local', 'threshold', 'maxiter']

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
            log.cite('bultinck2007', 'the use of Hirshfeld-I partitioning')

    def _proatom_change(self, number, old, new):
        '''Compute the difference between an old and a new spline

           **Arguments:**

           number
                Element number

           old, new
                The old and new splines
        '''
        # TODO: use the same criterion in CPart
        rtf = self.proatomdb.get_rtransform(number)
        return (4*np.pi)*dot_multi(
            rtf.get_radii()**2, # TODO: get routines are slow
            rtf.get_volume_elements(),
            SimpsonIntegrator1D().get_weights(rtf.npoint),
            (old.copy_y() - new.copy_y())**2 # TODO: copy is slow
        )

    def _get_proatom_fn(self, index, number, target_charge, first, grid):
        icharge = int(np.floor(target_charge))
        x = target_charge - icharge
        # check if icharge record should exist
        pseudo_pop = self.system.pseudo_numbers[index] - icharge
        if pseudo_pop == 1:
            return self._proatomdb.get_spline(number, {icharge: 1-x})
        elif pseudo_pop > 1:
            return self._proatomdb.get_spline(number, {icharge: 1-x, icharge+1: x})
        else:
            raise ValueError('Requesting a pro-atom with a negative (pseudo) population')

    @just_once
    def _init_at_weights(self):
        # Perform one general check in the beginning to keep things simple.
        if all(self.cache.has('at_weights', i) for i in xrange(self.system.natom)):
            return

        self.do_mol_dens()

        # Iterative Hirshfeld loop
        charges = np.zeros(self.system.natom)
        pseudo_populations = np.zeros(self.system.natom)
        counter = 0
        if log.do_medium:
            log('Iterative Hirshfeld partitioning loop')
            log.hline()
            log('Counter      Change   QTot Error')
            log.hline()

        while True:
            # Update current pro-atoms
            change = 0.0
            first_iter = True
            for i, grid in self.iter_grids():
                old_proatom_fn = self.cache.load('proatom_fn', i, default=None)
                number = self.system.numbers[i]
                charge = charges[i]
                proatom_fn = self._get_proatom_fn(i, number, charge, old_proatom_fn is None, grid)

                self.cache.dump('proatom_fn', i, proatom_fn)
                if old_proatom_fn is not None:
                    first_iter = False
                    change += self._proatom_change(number, old_proatom_fn, proatom_fn)

            # Enforce (single) update of pro-molecule in case of a global grid
            if not self.local:
                self.cache.invalidate('promol', self._molgrid.size)

            # Compute populations
            for i, grid in self.iter_grids():
                # Compute weight
                at_weights, new = self.cache.load('at_weights', i, alloc=grid.size)
                self._compute_at_weights(i, grid, at_weights)

                # Compute population
                dens = self.cache.load('mol_dens', i)
                pseudo_populations[i] = grid.integrate(at_weights, dens)
            charges = self.system.pseudo_numbers - pseudo_populations

            change = np.sqrt(change/self.system.natom)
            if log.do_medium:
                pop_error = charges.sum() - self.system.charge
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
        populations = pseudo_populations - (self.system.pseudo_numbers - self.system.numbers)
        self.cache.dump('pseudo_populations', pseudo_populations)
        self.cache.dump('populations', populations)
        self.cache.dump('charges', charges)
        self.cache.dump('niter', counter)
        self.cache.dump('change', change)

    def do_all(self):
        '''Computes all AIM properties and returns a corresponding list of keys'''
        names = HirshfeldDPart.do_all(self)
        return names + ['niter', 'change']
