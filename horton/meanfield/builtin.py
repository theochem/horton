# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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
'''Built-in energy terms'''


import numpy as np

from horton.meanfield.gridgroup import GridObservable
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = ['DiracExchange']


class DiracExchange(GridObservable):
    '''An implementation of the Dirac Exchange Functional'''
    def __init__(self, wfn, label='exchange_dirac', coeff=None):
        '''
           **Arguments:**

           wfn
                A Wavefunction object.

           **Optional arguments:**

           label
                A label for this observable.

           coeff
                The coefficient Cx in front of the Dirac exchange energy.
                It defaults to the uniform electron gas value, i.e.
                Cx = 3/4 (3/pi)^(1/3).
        '''
        self._wfn = wfn
        if coeff is None:
            self.coeff = 3.0 / 4.0 * (3.0 / np.pi) ** (1.0 / 3.0)
        else:
            self.coeff = coeff
        self.derived_coeff = -self.coeff * (4.0 / 3.0) * 2 ** (1.0 / 3.0)
        GridObservable.__init__(self, label)

    def compute_energy(self, cache, grid):
        self._update_pot(cache, grid)

        def helper(select):
            pot = cache['pot_exchange_dirac_%s' % select]
            rho = cache['rho_%s' % select]
            return grid.integrate(pot, rho)

        energy = helper('alpha')
        if isinstance(self._wfn, RestrictedWFN):
            energy *= 2
        else:
            energy += helper('beta')
        energy *= 3.0 / 4.0
        return energy

    def _update_pot(self, cache, grid):
        '''Recompute the Exchange potential(s) if invalid'''
        def helper(select):
            # update grid stuff
            rho = cache['rho_%s' % select]
            pot, new = cache.load('pot_exchange_dirac_%s' % select, alloc=grid.size)
            if new:
                pot[:] = self.derived_coeff * rho ** (1.0 / 3.0)

        helper('alpha')
        if isinstance(self._wfn, UnrestrictedWFN):
            helper('beta')

    def compute_pot(self, cache, grid):
        '''Recompute the Exchange potential(s) if invalid'''
        self._update_pot(cache, grid)
        cache['dpot_total_alpha'] += cache['pot_exchange_dirac_alpha']
        if isinstance(self._wfn, UnrestrictedWFN):
            cache['dpot_total_beta'] += cache['pot_exchange_dirac_beta']
