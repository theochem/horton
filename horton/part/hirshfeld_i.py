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
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import dot_multi
from horton.log import log
from horton.part.hirshfeld import HirshfeldWPart, HirshfeldCPart


__all__ = ['HirshfeldIWPart', 'HirshfeldICPart']


class HirshfeldIMixin(object):
    name = 'hi'
    options = ['threshold', 'maxiter']

    def __init__(self, threshold=1e-6, maxiter=500):
        self._threshold = threshold
        self._maxiter = maxiter

    def _init_log_scheme(self):
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld-I'),
                ('Convergence threshold', '%.1e' % self._threshold),
                ('Maximum iterations', self._maxiter),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('bultinck2007', 'the use of Hirshfeld-I partitioning')

    def get_interpolation_info(self, i, charges=None):
        if charges is None:
            charges = self.cache.load('charges')
        target_charge = charges[i]
        icharge = int(np.floor(target_charge))
        x = target_charge - icharge
        return icharge, x

    def get_proatom_rho(self, index, charges=None):
        icharge, x = self.get_interpolation_info(index, charges)
        # check if icharge record should exist
        pseudo_pop = self.system.pseudo_numbers[index] - icharge
        number = self.system.numbers[index]
        if pseudo_pop == 1:
            return self.proatomdb.get_rho(number, {icharge: 1-x})
        elif pseudo_pop > 1:
            return self.proatomdb.get_rho(number, {icharge: 1-x, icharge+1: x})
        else:
            raise ValueError('Requesting a pro-atom with a negative (pseudo) population')

    def compute_change(self, propars1, propars2):
        '''Compute the difference between an old and a new proatoms'''
        msd = 0.0 # mean-square deviations
        for i in xrange(self.system.natom):
            number = self.system.numbers[i]
            weights = self.proatomdb.get_radial_weights(number)
            rho1 = self.get_proatom_rho(i, propars1)
            rho2 = self.get_proatom_rho(i, propars2)
            delta = rho1 - rho2
            msd +=  dot_multi(weights, delta, delta)
        return np.sqrt(msd)

    def _init_propars(self):
        self.history_propars = []
        self.history_charges = []
        charges = self.cache.load('charges', alloc=self._system.natom)[0]
        self.cache.dump('propars', charges)
        return charges

    def _update_propars(self):
        # Keep track of history
        self.history_propars.append(self.cache.load('propars').copy())

        # Enforce (single) update of pro-molecule in case of a global grid
        self.cache.invalidate('promoldens')

        for index in xrange(self.natom):
            # Update proatom
            self._update_propars_atom(index)

        # Keep track of history
        self.history_charges.append(self.cache.load('charges').copy())

    def _update_propars_atom(self, index):
        # Compute population
        self.cache.invalidate('at_weights', index)
        pseudo_population = self.compute_pseudo_population(index)

        # Store charge
        charges = self.cache.load('charges')
        charges[index] = self.system.pseudo_numbers[index] - pseudo_population

    def _finalize_propars(self):
        charges = self._cache.load('charges')
        self.cache.dump('history_propars', np.array(self.history_propars))
        self.cache.dump('history_charges', np.array(self.history_charges))
        self.cache.dump('populations', self.system.numbers - charges)
        self.cache.dump('pseudo_populations', self.system.pseudo_numbers - charges)

    @just_once
    def _init_partitioning(self):
        # Perform one general check in the beginning to keep things simple.
        if all(self.cache.has('at_weights', i) for i in xrange(self.system.natom)):
            return
        # Need to compute density
        self.do_moldens()

        propars = self._init_propars()
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()

        counter = 0
        change = 1e100

        while True:
            counter += 1

            # Update the parameters that determine the pro-atoms.
            old_propars = propars.copy()
            self._update_propars()

            # Check for convergence
            change = self.compute_change(propars, old_propars)
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < self._threshold or counter >= self._maxiter:
                break

        if log.medium:
            log.hline()

        self._finalize_propars()
        self.cache.dump('niter', counter)
        self.cache.dump('change', change)


class HirshfeldIWPart(HirshfeldIMixin, HirshfeldWPart):
    def __init__(self, system, grid, proatomdb, local=True, threshold=1e-6, maxiter=500):
        '''
           **Optional arguments:** (that are not present in the base class)

           threshold
                The procedure is considered to be converged when the maximum
                change of the charges between two iterations drops below this
                threshold.

           maxiter
                The maximum number of iterations. If no convergence is reached
                in the end, no warning is given.
        '''
        HirshfeldIMixin.__init__(self, threshold, maxiter)
        HirshfeldWPart.__init__(self, system, grid, proatomdb, local)

    def do_all(self):
        '''Computes all AIM properties and returns a corresponding list of keys'''
        names = HirshfeldWPart.do_all(self)
        return names + ['niter', 'change', 'history_charges', 'history_propars']


class HirshfeldICPart(HirshfeldIMixin, HirshfeldCPart):
    def __init__(self, system, grid, local, moldens, proatomdb, wcor_numbers, wcor_rcut_max=2.0, wcor_rcond=0.1, threshold=1e-6, maxiter=500):
        '''
           **Optional arguments:** (that are not present in the base class)

           threshold
                The procedure is considered to be converged when the maximum
                change of the charges between two iterations drops below this
                threshold.

           maxiter
                The maximum number of iterations. If no convergence is reached
                in the end, no warning is given.
        '''
        HirshfeldIMixin.__init__(self, threshold, maxiter)
        HirshfeldCPart.__init__(self, system, grid, local, moldens, proatomdb, wcor_numbers, wcor_rcut_max, wcor_rcond)

    def do_all(self):
        names = HirshfeldCPart.do_all(self)
        return names + ['niter', 'change', 'history_charges', 'history_propars']
