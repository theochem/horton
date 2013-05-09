# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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


import numpy as np

from horton.cache import just_once
from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import dot_multi
from horton.log import log
from horton.part.hirshfeld import HirshfeldWPart, HirshfeldCPart


__all__ = ['HirshfeldIWPart', 'HirshfeldICPart']


class HirshfeldIMixin(object):
    name = 'hi'
    options = ['threshold', 'maxiter', 'greedy']

    def __init__(self, threshold=1e-6, maxiter=500, greedy=False):
        self._threshold = threshold
        self._maxiter = maxiter
        self._greedy = greedy

    def _init_log_scheme(self):
        if log.do_medium:
            log.deflist([
                ('Scheme', 'Hirshfeld-I'),
                ('Convergence threshold', '%.1e' % self._threshold),
                ('Maximum iterations', self._maxiter),
                ('Proatomic DB',  self._proatomdb),
            ])
            log.cite('bultinck2007', 'the use of Hirshfeld-I partitioning')

    def get_memory_estimates(self):
        if self._greedy:
            return [('Isolated atoms', np.ones(self.natom)*3, 0),] # This is a conservative estimate.
        else:
            return []

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
        if pseudo_pop == 1 or x == 0.0:
            return self.proatomdb.get_rho(number, {icharge: 1-x})
        elif pseudo_pop > 1:
            return self.proatomdb.get_rho(number, {icharge: 1-x, icharge+1: x})
        elif psuedo_pop <= 0:
            raise ValueError('Requesting a pro-atom with a negative (pseudo) population')

    def compute_change(self, propars1, propars2):
        '''Compute the difference between an old and a new proatoms'''
        msd = 0.0 # mean-square deviation
        for i in xrange(self.system.natom):
            number = self.system.numbers[i]
            rgrid = self.proatomdb.get_rgrid(number)
            rho1 = self.get_proatom_rho(i, propars1)
            rho2 = self.get_proatom_rho(i, propars2)
            delta = rho1 - rho2
            msd +=  rgrid.integrate(delta, delta)
        return np.sqrt(msd)

    def get_somefn(self, index, spline, key, label, grid=None):
        if grid is None:
            grid = self.get_grid(index)
        key = key + (index, id(grid))
        if self._greedy:
            result, new = self.cache.load(*key, alloc=grid.shape)
        else:
            result = grid.zeros()
            new = True
        if new:
            self.eval_spline(index, spline, result, grid, label)
        return result

    def get_isolated(self, index, charge, grid=None):
        number = self.system.numbers[index]
        spline = self.proatomdb.get_spline(number, charge)
        return self.get_somefn(index, spline, ('isolated', charge), 'isolated q=%+i' % charge, grid)

    def eval_proatom(self, index, output, grid=None):
        # Greedy version of eval_proatom
        icharge, x = self.get_interpolation_info(index)
        output[:] = self.get_isolated(index, icharge, grid)
        output *= 1-x
        pseudo_pop = self.system.pseudo_numbers[index] - icharge
        if pseudo_pop > 1 and x != 0.0:
            output += self.get_isolated(index, icharge+1, grid)*x
        elif pseudo_pop <= 0:
            raise ValueError('Requesting a pro-atom with a negative (pseudo) population')
        output += 1e-100

    def _init_propars(self):
        self.history_propars = []
        self.history_charges = []
        charges = self.cache.load('charges', alloc=self._system.natom)[0]
        self.cache.dump('propars', charges)
        return charges

    def _update_propars(self):
        # Keep track of history
        self.history_propars.append(self.cache.load('propars').copy())

        # Update the partitioning based on the latest proatoms
        self.update_at_weights()

        # Update the proatoms
        for index in xrange(self.natom):
            self._update_propars_atom(index)

        # Keep track of history
        self.history_charges.append(self.cache.load('charges').copy())

    def _update_propars_atom(self, index):
        # Compute population
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
    def do_partitioning(self):
        # Perform one general check in the beginning to avoid recomputation
        new = any(('at_weights', i) not in self.cache for i in xrange(self.system.natom))
        new |= 'niter' not in self.cache
        new |= 'change'not in self.cache
        if new:
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
    do_partitioning.names = ['niter', 'change', 'history_charges', 'history_propars']


class HirshfeldIWPart(HirshfeldIMixin, HirshfeldWPart):
    def __init__(self, system, grid, proatomdb, local=True, threshold=1e-6, maxiter=500, greedy=False):
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
        HirshfeldIMixin.__init__(self, threshold, maxiter, greedy)
        HirshfeldWPart.__init__(self, system, grid, proatomdb, local)

    def get_memory_estimates(self):
        return (
            HirshfeldWPart.get_memory_estimates(self) +
            HirshfeldIMixin.get_memory_estimates(self)
        )

    def eval_proatom(self, index, output, grid=None):
        if self._greedy:
            HirshfeldIMixin.eval_proatom(self, index, output, grid)
        else:
            HirshfeldWPart.eval_proatom(self, index, output, grid)


class HirshfeldICPart(HirshfeldIMixin, HirshfeldCPart):
    def __init__(self, system, grid, local, moldens, proatomdb, wcor_numbers, wcor_rcut_max=2.0, wcor_rcond=0.1, threshold=1e-6, maxiter=500, greedy=False):
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
        HirshfeldIMixin.__init__(self, threshold, maxiter, greedy)
        HirshfeldCPart.__init__(self, system, grid, local, moldens, proatomdb, wcor_numbers, wcor_rcut_max, wcor_rcond)

    def get_memory_estimates(self):
        return (
            HirshfeldCPart.get_memory_estimates(self) +
            HirshfeldIMixin.get_memory_estimates(self)
        )

    def eval_proatom(self, index, output, grid=None):
        if self._greedy:
            HirshfeldIMixin.eval_proatom(self, index, output, grid)
        else:
            HirshfeldCPart.eval_proatom(self, index, output, grid)
