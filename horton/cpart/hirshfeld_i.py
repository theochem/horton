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

from horton.cpart.hirshfeld import HirshfeldCPart
from horton.log import log


__all__ = [
    'HirshfeldICPart',
]


# TODO: reduce duplicate code

class HirshfeldICPart(HirshfeldCPart):
    name = 'hi'
    options = ['smooth', 'max_iter', 'threshold']

    def __init__(self, system, ui_grid, moldens, proatomdb, store, smooth=False, max_iter=100, threshold=1e-4):
        '''
           **Optional arguments:** (those not present in the base class)

           max_iter
                The maximum number of iterations. If no convergence is reached
                in the end, no warning is given.

           threshold
                The procedure is considered to be converged when the maximum
                change of the charges between two iterations drops below this
                threshold.
        '''
        self._max_iter = max_iter
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

    def _init_propars(self):
        self.history = []
        return self._cache.load('charges', alloc=self._system.natom)[0]

    def _update_propars(self, charges, work):
        self.history.append(charges.copy())
        for index in xrange(self._system.natom):
            pseudo_population = self.compute_pseudo_population(index, work)
            charges[index] = self.system.pseudo_numbers[index] - pseudo_population

    def _finalize_propars(self, charges):
        self.history.append(charges.copy())
        self._cache.dump('history_charges', np.array(self.history))
        self._cache.dump('populations', self.system.numbers - charges)

    def _init_partitioning(self):
        propars = self._init_propars()
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()

        counter = 0
        change = 1e100
        work = self._ui_grid.zeros()

        while True:
            counter += 1

            # Update the pro-molecule density
            self._update_promolecule(work)

            # Compute the atomic weight functions if this is useful. This is merely
            # a matter of efficiency.
            self._store_at_weights(work)

            # Update the parameters that determine the pro-atoms.
            old_propars = propars.copy()
            self._update_propars(propars, work)

            # Check for convergence
            change = abs(propars - old_propars).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < self._threshold or counter >= self._max_iter:
                break

        if log.medium:
            log.hline()

        self._finalize_propars(propars)
        self._cache.dump('niter', counter)
        self._cache.dump('change', change)

    def get_interpolation_info(self, i):
        target_charge = self._cache.load('charges')[i]
        icharge = int(np.floor(target_charge))
        x = target_charge - icharge
        return icharge, x

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

    def get_proatom_spline(self, index):
        icharge, x = self.get_interpolation_info(index)
        number = self._system.numbers[index]
        return self._proatomdb.get_spline(number, {icharge: 1-x, icharge+1: x})

    def do_all(self):
        names = HirshfeldCPart.do_all(self)
        return names + ['niter', 'change', 'history_charges']
