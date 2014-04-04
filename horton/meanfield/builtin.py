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


__all__ = ['RDiracExchange', 'UDiracExchange']


class DiracExchange(GridObservable):
    '''Common code for the Dirac Exchange Functional implementations'''
    def __init__(self, label='x_dirac', coeff=None):
        '''
           **Optional arguments:**

           label
                A label for this observable.

           coeff
                The coefficient Cx in front of the Dirac exchange energy.
                It defaults to the uniform electron gas value, i.e.
                Cx = 3/4 (3/pi)^(1/3).
        '''
        if coeff is None:
            self.coeff = 3.0 / 4.0 * (3.0 / np.pi) ** (1.0 / 3.0)
        else:
            self.coeff = coeff
        self.derived_coeff = -self.coeff * (4.0 / 3.0) * 2 ** (1.0 / 3.0)
        GridObservable.__init__(self, label)

    def _update_pot(self, cache, grid, select):
        '''Recompute an Exchange potential if invalid

           cache
                A cache instance in which the potential is stored.

           grid
                The integration grid.

           select
                'alpha' or 'beta'
        '''
        rho = cache['rho_%s' % select]
        pot, new = cache.load('pot_x_dirac_%s' % select, alloc=grid.size)
        if new:
            pot[:] = self.derived_coeff * rho ** (1.0 / 3.0)
        return pot


class RDiracExchange(DiracExchange):
    '''The Dirac Exchange Functional for restricted wavefunctions'''

    def compute(self, cache, grid):
        '''See ``GridObservable.compute``'''
        pot = self._update_pot(cache, grid, 'alpha')
        rho = cache['rho_alpha']
        return (3.0 / 2.0) * grid.integrate(pot, rho)

    def add_pot(self, cache, grid, dpot_alpha):
        '''See ``GridObservable.add_pot``'''
        dpot_alpha += self._update_pot(cache, grid, 'alpha')


class UDiracExchange(DiracExchange):
    '''The Dirac Exchange Functional for unrestricted wavefunctions'''

    def compute(self, cache, grid):
        '''See ``GridObservable.compute``'''
        pot_alpha = self._update_pot(cache, grid, 'alpha')
        pot_beta = self._update_pot(cache, grid, 'beta')
        rho_alpha = cache['rho_alpha']
        rho_beta = cache['rho_beta']
        return (3.0 / 4.0) * (grid.integrate(pot_alpha, rho_alpha) +
                              grid.integrate(pot_beta, rho_beta))

    def add_pot(self, cache, grid, dpot_alpha, dpot_beta):
        '''See ``GridObservable.add_pot``'''
        dpot_alpha += self._update_pot(cache, grid, 'alpha')
        dpot_beta += self._update_pot(cache, grid, 'beta')
