# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
'''Built-in energy terms'''


import numpy as np

from horton.meanfield.gridgroup import GridObservable
from horton.grid.molgrid import BeckeMolGrid
from horton.grid.poisson import solve_poisson_becke
from horton.utils import doc_inherit


__all__ = ['RBeckeHartree', 'UBeckeHartree', 'RDiracExchange', 'UDiracExchange']


class BeckeHartree(GridObservable):
    def __init__(self, lmax, label='hartree_becke'):
        self.lmax = lmax
        GridObservable.__init__(self, label)

    def _update_pot(self, cache, grid):
        '''Recompute an Hartree potential if invalid

           cache
                A cache instance in which the potential is stored.

           grid
                The integration grid.
        '''
        # This only works under a few circumstances
        if not isinstance(grid, BeckeMolGrid):
            raise TypeError('The BeckeHatree term only works for Becke-Lebedev molecular integration grids')
        if grid.mode != 'keep':
            raise TypeError('The mode option of the molecular grid must be \'keep\'.')

        pot, new = cache.load('pot_%s' % self.label, alloc=grid.size)
        if new:
            rho = cache['rho_full']
            # Construct spherical decompositions of atomic densities, derive
            # hartree potentials and evaluate
            begin = 0
            pot[:] = 0
            for atgrid in grid.subgrids:
                end = begin + atgrid.size
                becke_weights = grid.becke_weights[begin:end]
                density_decomposition = atgrid.get_spherical_decomposition(rho[begin:end], becke_weights, lmax=self.lmax)
                hartree_decomposition = solve_poisson_becke(density_decomposition)
                grid.eval_decomposition(hartree_decomposition, atgrid.center, pot)
                begin = end
        return pot

    @doc_inherit(GridObservable)
    def compute_energy(self, cache, grid):
        pot = self._update_pot(cache, grid)
        rho = cache['rho_full']
        return 0.5*grid.integrate(pot, rho)


class RBeckeHartree(BeckeHartree):
    @doc_inherit(BeckeHartree)
    def add_pot(self, cache, grid, dpot_alpha):
        pot = self._update_pot(cache, grid)
        dpot_alpha += pot


class UBeckeHartree(BeckeHartree):
    @doc_inherit(BeckeHartree)
    def add_pot(self, cache, grid, dpot_alpha, dpot_beta):
        pot = self._update_pot(cache, grid)
        dpot_alpha += pot
        dpot_beta += pot


class DiracExchange(GridObservable):
    '''Common code for the Dirac Exchange Functional implementations'''
    def __init__(self, label='x_dirac', coeff=None):
        r'''
           **Optional arguments:**

           label
                A label for this observable.

           coeff
                The coefficient Cx in front of the Dirac exchange energy.
                It defaults to the uniform electron gas value, i.e.
                :math:`C_x = \frac{3}{4} \left(\frac{3}{\pi}\right)^{1/3}`.
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

    @doc_inherit(GridObservable)
    def compute_energy(self, cache, grid):
        pot = self._update_pot(cache, grid, 'alpha')
        rho = cache['rho_alpha']
        return (3.0 / 2.0) * grid.integrate(pot, rho)

    @doc_inherit(GridObservable)
    def add_pot(self, cache, grid, dpot_alpha):
        dpot_alpha += self._update_pot(cache, grid, 'alpha')


class UDiracExchange(DiracExchange):
    '''The Dirac Exchange Functional for unrestricted wavefunctions'''

    @doc_inherit(GridObservable)
    def compute_energy(self, cache, grid):
        pot_alpha = self._update_pot(cache, grid, 'alpha')
        pot_beta = self._update_pot(cache, grid, 'beta')
        rho_alpha = cache['rho_alpha']
        rho_beta = cache['rho_beta']
        return (3.0 / 4.0) * (grid.integrate(pot_alpha, rho_alpha) +
                              grid.integrate(pot_beta, rho_beta))

    @doc_inherit(GridObservable)
    def add_pot(self, cache, grid, dpot_alpha, dpot_beta):
        dpot_alpha += self._update_pot(cache, grid, 'alpha')
        dpot_beta += self._update_pot(cache, grid, 'beta')
