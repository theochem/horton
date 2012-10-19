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

from horton.hamiltonian.core import HamiltonianTerm
from horton.hamiltonian.cext import LibXCWrapper


__all__ = ['LibXCTerm']


class LibXCTerm(HamiltonianTerm):
    '''Any LDA functional from LibXC'''

    require_grid = True
    def __init__(self, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC
        '''
        self._name = name
        self._libxc_wrapper = LibXCWrapper(name)

    def _update_operator(self):
        if self.system.wfn.closed_shell:
            # In the closed-shell case, libxc expects the total density as input
            # and returns the potential for the alpha electrons.
            pot, new = self.cache.load('pot_libxc_%s_alpha' % self._name, alloc=self.grid.size)
            if new:
                rho = self.update_rho('full')
                self._libxc_wrapper.compute_vxc_unpol(rho, pot)

            # TODO: in the Hamiltonian class, all the grids should be added
            # and converted only once to a one-body operator
            # Create/update the one-body operator based on the potential on the grid.
            operator, new = self.cache.load('op_libxc_%s_alpha' % self._name, alloc=(self.system.lf, 'one_body'))
            if new:
                self.system.compute_grid_one_body(self.grid.points, self.grid.weights, pot, operator)
        else:
            # In case of spin-polarized computations, alpha and beta densities
            # go in and the alpha and beta potentials come out of  it.
            pot_alpha, new_alpha = self.cache.load('pot_libxc_%s_alpha' % self._name, alloc=self.grid.size)
            pot_beta, new_beta = self.cache.load('pot_libxc_%s_beta' % self._name, alloc=self.grid.size)
            if new_alpha or new_beta:
                # TODO: rho_both has to be cashed
                rho_both = np.zeros((len(pot_alpha), 2), float)
                rho_both[:,0] = self.update_rho('alpha')
                rho_both[:,1] = self.update_rho('beta')
                # TODO: pot_both has to be cashed
                pot_both = np.zeros((len(pot_alpha), 2), float)
                self._libxc_wrapper.compute_vxc_pol(rho_both, pot_both)
                pot_alpha[:] = pot_both[:,0]
                pot_beta[:] = pot_both[:,1]

            # TODO: in the Hamiltonian class, all the grids should be added
            # and converted only once to a one-body operator
            # Create/update the one-body operator based on the potential on the grid.
            operator_alpha, new_alpha = self.cache.load('op_libxc_%s_alpha' % self._name, alloc=(self.system.lf, 'one_body'))
            if new_alpha:
                self.system.compute_grid_one_body(self.grid.points, self.grid.weights, pot_alpha, operator_alpha)

            operator_beta, new_beta = self.cache.load('op_libxc_%s_beta' % self._name, alloc=(self.system.lf, 'one_body'))
            if new_beta:
                self.system.compute_grid_one_body(self.grid.points, self.grid.weights, pot_beta, operator_beta)

    def compute_energy(self):
        if self.system.wfn.closed_shell:
            # In the unpolarized case, libxc expects the total density as input
            # and returns the energy density per particle.
            rho = self.update_rho('full')
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                self._libxc_wrapper.compute_exc_unpol(rho, edens)
            energy = self.grid.integrate(edens, rho)
            self.store_energy('libxc_%s' % self._name, energy)
            return energy
        else:
            # In case of spin-polarized computations, alpha and beta densities
            # go in and the 'total' energy density comes out.
            rho_alpha = self.update_rho('alpha')
            rho_beta = self.update_rho('beta')
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                # TODO: rho_both has to be cashed
                rho_both = np.zeros((len(rho_alpha), 2), float)
                rho_both[:,0] = rho_alpha
                rho_both[:,1] = rho_beta
                self._libxc_wrapper.compute_exc_pol(rho_both, edens)

            energy = self.grid.integrate(edens, rho_alpha) + self.grid.integrate(edens, rho_beta)
            self.store_energy('libxc_%s' % self._name, energy)
            return energy

    def add_fock_matrix(self, fock_alpha, fock_beta):
        self._update_operator()
        fock_alpha.iadd(self.cache.load('op_libxc_%s_alpha' % self._name))
        if not self.system.wfn.closed_shell:
            fock_beta.iadd(self.cache.load('op_libxc_%s_beta' % self._name, ))
