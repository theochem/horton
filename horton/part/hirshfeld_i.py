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
from horton.part.base import DPart
from horton.part.hirshfeld import HirshfeldDPart, HirshfeldCPart


__all__ = ['HirshfeldIDPart', 'HirshfeldICPart']


# TODO: isolate common code in mixin class
# TODO: proofread and add tests for pseudo densities


class HirshfeldIMixin(object):
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
        return self.cache.load('charges', alloc=self._system.natom)[0]

    def _finalize_propars(self, propars):
        charges = self._cache.load('charges')
        self.cache.dump('history_propars', np.array(self.history_propars))
        self.cache.dump('history_charges', np.array(self.history_charges))
        self.cache.dump('populations', self.system.numbers - charges)
        self.cache.dump('pseudo_populations', self.system.pseudo_numbers - charges)

    @just_once
    def _init_partitioning(self):
        # TODO: do not pass around propars array
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
            self._update_propars(propars)

            # Check for convergence
            change = self.compute_change(propars, old_propars)
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < self._threshold or counter >= self._maxiter:
                break

        if log.medium:
            log.hline()

        self._finalize_propars(propars)
        self.cache.dump('niter', counter)
        self.cache.dump('change', change)


class HirshfeldIDPart(HirshfeldIMixin, HirshfeldDPart):
    '''Iterative Hirshfeld partitioning'''

    name = 'hi'
    options = ['local', 'threshold', 'maxiter']

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

    def _update_propars(self, propars):
        # Keep track of history
        self.history_propars.append(propars.copy())

        # Enforce (single) update of pro-molecule in case of a global grid
        if not self.local:
            self.cache.invalidate('promol', self._molgrid.size)

        for index, grid in self.iter_grids():# TODO: Get rid of this construct
            # Update proatom
            self._update_propars_atom(index, grid, propars)

        # Keep track of history
        self.history_charges.append(self._cache.load('charges').copy())

    def _update_propars_atom(self, index, grid, propars):
        # Compute weight
        at_weights, new = self.cache.load('at_weights', index, alloc=grid.size)
        self.compute_at_weights(index, grid, at_weights)

        # Compute population
        dens = self.cache.load('mol_dens', index)
        charges = self._cache.load('charges')
        pseudo_population = grid.integrate(at_weights, dens)
        charges[index] = self.system.pseudo_numbers[index] - pseudo_population

    @just_once
    def _init_partitioning(self):
        # Perform one general check in the beginning to keep things simple.
        if all(self.cache.has('at_weights', i) for i in xrange(self.system.natom)):
            return
        # Need to compute density
        self.do_mol_dens()
        # Run the generic code
        HirshfeldIMixin._init_partitioning(self)


    def do_all(self):
        '''Computes all AIM properties and returns a corresponding list of keys'''
        names = HirshfeldDPart.do_all(self)
        return names + ['niter', 'change', 'history_charges', 'history_propars']


class HirshfeldICPart(HirshfeldIMixin, HirshfeldCPart):
    '''Iterative Hirshfeld partitioning'''

    name = 'hi'
    options = ['smooth', 'maxiter', 'threshold']

    def __init__(self, system, ui_grid, moldens, proatomdb, store, smooth=False, maxiter=100, threshold=1e-4):
        '''
           **Optional arguments:** (those not present in the base class)

           maxiter
                The maximum number of iterations. If no convergence is reached
                in the end, no warning is given.

           threshold
                The procedure is considered to be converged when the maximum
                change of the charges between two iterations drops below this
                threshold.
        '''
        self._maxiter = maxiter
        self._threshold = threshold
        HirshfeldCPart.__init__(self, system, ui_grid, moldens, proatomdb, store, smooth)

    def _get_isolated_atom(self, i, charge, output):
        key = ('isolated_atom', i, charge)
        if key in self._store:
            self._store.load(output, *key)
        else:
            number = self._system.numbers[i]
            spline = self._proatomdb.get_spline(number, charge)
            self.compute_spline(i, spline, output, 'n=%i q=%+i' % (number, charge))
            self._store.dump(output, *key)

    def compute_proatom(self, i, output, window=None):
        if self._store.fake or window is not None:
            HirshfeldCPart.compute_proatom(self, i, output, window)
        else:
            # Construct the pro-atom
            icharge, x = self.get_interpolation_info(i)
            self._get_isolated_atom(i, icharge, output)
            if x != 1:
                output *= 1-x
                ipseudo_pop = self.system.pseudo_numbers[i] - icharge
                if ipseudo_pop > 1:
                    tmp = self._ui_grid.zeros()
                    self._get_isolated_atom(i, icharge+1, tmp)
                    tmp *= x
                    output += tmp
            output += 1e-100 # avoid division by zero

    def _update_propars_atom(self, index, charges, work):
        pseudo_population = self.compute_pseudo_population(index, work)
        charges[index] = self.system.pseudo_numbers[index] - pseudo_population

    def _update_propars(self, propars):
        self.history_propars.append(propars.copy())

        # TODO: avoid periodic reallocation of work array
        work = self._ui_grid.zeros()

        # Update the pro-molecule density
        self._update_promolecule(work)

        # Compute the atomic weight functions if this is useful. This is merely
        # a matter of efficiency.
        self._store_at_weights(work)

        # Update the pro-atom parameters.
        for index in xrange(self._system.natom):
            self._update_propars_atom(index, propars, work)

        self.history_charges.append(self.cache.load('charges').copy())

    def do_all(self):
        names = HirshfeldCPart.do_all(self)
        return names + ['niter', 'change', 'history_charges', 'history_propars']
