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


import numpy as np

from horton.log import log
from horton.meanfield.observable import Observable
from horton.meanfield.cext import LibXCWrapper


__all__ = ['LibXCLDA', 'LibXCGGA', 'LibXCHybridGGA']



class LibXCEnergy(Observable):
    def __init__(self, name):
        self._name = name
        self._libxc_wrapper = LibXCWrapper(name)
        log.cite('marques2012', 'using LibXC, the library of exchange and correlation functionals')

    def _update_operator(self):
        raise NotImplementedError

    def add_fock_matrix(self, fock_alpha, fock_beta):
        # TODO: move this above compute, also in all other classes
        self._update_operator()
        fock_alpha.iadd(self.cache.load('op_libxc_%s_alpha' % self._name))
        if not self.system.wfn.closed_shell:
            fock_beta.iadd(self.cache.load('op_libxc_%s_beta' % self._name, ))


class LibXCLDA(LibXCEnergy):
    '''Any LDA functional from LibXC'''

    require_grid = True
    def __init__(self, name):
        '''
           **Arguments:*

           name
                The name of the functional in LibXC, without the 'lda_' prefix.
        '''
        LibXCEnergy.__init__(self, 'lda_' + name.lower())

    def _update_operator(self):
        if self.system.wfn.closed_shell:
            # In the closed-shell case, libxc expects the total density as input
            # and returns the potential for the alpha electrons.
            pot, new = self.cache.load('pot_libxc_%s_alpha' % self._name, alloc=self.grid.size)
            if new:
                rho = self.update_rho('full')
                self._libxc_wrapper.compute_lda_vxc_unpol(rho, pot)

            # TODO: in the Hamiltonian class, all the grids should be added
            # and converted only once to a fock operator
            # Create/update the one-body operator based on the potential on the grid.
            operator, new = self.cache.load('op_libxc_%s_alpha' % self._name, alloc=(self.system.lf, 'one_body'))
            if new:
                self.system.compute_grid_density_fock(self.grid.points, self.grid.weights, pot, operator)
        else:
            # In case of spin-polarized computations, alpha and beta densities
            # go in and the alpha and beta potentials come out of  it.
            pot_both, new = self.cache.load('pot_libxc_%s_both' % self._name, alloc=(self.grid.size, 2))
            if new:
                rho_both = self.update_rho('both')
                self._libxc_wrapper.compute_lda_vxc_pol(rho_both, pot_both)

            # TODO: in the Hamiltonian class, all the grids should be added
            # and converted only once to a fock operator
            # Create/update the one-body operator based on the potential on the grid.
            operator_alpha, new_alpha = self.cache.load('op_libxc_%s_alpha' % self._name, alloc=(self.system.lf, 'one_body'))
            if new_alpha:
                self.system.compute_grid_density_fock(self.grid.points, self.grid.weights, pot_both[:,0], operator_alpha)

            operator_beta, new_beta = self.cache.load('op_libxc_%s_beta' % self._name, alloc=(self.system.lf, 'one_body'))
            if new_beta:
                self.system.compute_grid_density_fock(self.grid.points, self.grid.weights, pot_both[:,1], operator_beta)

    def compute(self):
        if self.system.wfn.closed_shell:
            # In the unpolarized case, libxc expects the total density as input
            # and returns the energy density per particle.
            rho = self.update_rho('full')
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                self._libxc_wrapper.compute_lda_exc_unpol(rho, edens)
            energy = self.grid.integrate(edens, rho)
            self.store_energy('libxc_%s' % self._name, energy)
        else:
            # In case of spin-polarized computations, alpha and beta densities
            # go in and the 'total' energy density comes out.
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                rho_both = self.update_rho('both')
                self._libxc_wrapper.compute_lda_exc_pol(rho_both, edens)

            rho = self.update_rho('full')
            energy = self.grid.integrate(edens, rho)
            self.store_energy('libxc_%s' % self._name, energy)
        return energy


class LibXCGGA(LibXCEnergy):
    '''Any GGA functional from LibXC'''
    def __init__(self, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the 'gga_' prefix.
        '''
        LibXCEnergy.__init__(self, 'gga_' + name.lower())

    def _update_operator(self):
        if self.system.wfn.closed_shell:
            dpot, newd = self.cache.load('dpot_libxc_%s_alpha' % self._name, alloc=self.grid.size)
            spot, newt = self.cache.load('spot_libxc_%s_alpha' % self._name, alloc=self.grid.size)
            if newd or newt:
                rho = self.update_rho('full')
                sigma = self.update_sigma('full')
                self._libxc_wrapper.compute_gga_vxc_unpol(rho, sigma, dpot, spot)

            gpot, new = self.cache.load('gpot_libxc_%s_alpha' % self._name, alloc=(self.grid.size,3))
            if new:
                grad_rho = self.update_grad_rho('full')
                np.multiply(grad_rho, spot.reshape(-1,1), out=gpot)
                gpot *= 2


            # TODO: in the Hamiltonian class, all the grids should be added
            operator, new = self.cache.load('op_libxc_%s_alpha' % self._name, alloc=(self.system.lf, 'one_body'))
            if new:
                self.system.compute_grid_density_fock(self.grid.points, self.grid.weights, dpot, operator)
                self.system.compute_grid_gradient_fock(self.grid.points, self.grid.weights, gpot, operator)
        else:
            dpot_both, newd = self.cache.load('dpot_libxc_%s_both' % self._name, alloc=(self.grid.size, 2))
            spot_all, newt = self.cache.load('spot_libxc_%s_all' % self._name, alloc=(self.grid.size, 3))
            if newd or newt:
                rho_both = self.update_rho('both')
                sigma_all = self.update_sigma('all')
                self._libxc_wrapper.compute_gga_vxc_pol(rho_both, sigma_all, dpot_both, spot_all)

            gpot_alpha, new = self.cache.load('gpot_libxc_%s_alpha' % self._name, alloc=(self.grid.size,3))
            if new:
                # TODO: make more efficient
                gpot_alpha[:] = (2*spot_all[:,0].reshape(-1,1))*self.update_grad_rho('alpha')
                gpot_alpha[:] += (spot_all[:,1].reshape(-1,1))*self.update_grad_rho('beta')

            gpot_beta, new = self.cache.load('gpot_libxc_%s_beta' % self._name, alloc=(self.grid.size,3))
            if new:
                # TODO: make more efficient
                gpot_beta[:] = (2*spot_all[:,2].reshape(-1,1))*self.update_grad_rho('beta')
                gpot_beta[:] += (spot_all[:,1].reshape(-1,1))*self.update_grad_rho('alpha')


            # TODO: in the Hamiltonian class, all the grids should be added
            operator_alpha, new = self.cache.load('op_libxc_%s_alpha' % self._name, alloc=(self.system.lf, 'one_body'))
            if new:
                self.system.compute_grid_density_fock(self.grid.points, self.grid.weights, dpot_both[:,0], operator_alpha)
                self.system.compute_grid_gradient_fock(self.grid.points, self.grid.weights, gpot_alpha, operator_alpha)

            operator_beta, new = self.cache.load('op_libxc_%s_beta' % self._name, alloc=(self.system.lf, 'one_body'))
            if new:
                self.system.compute_grid_density_fock(self.grid.points, self.grid.weights, dpot_both[:,1], operator_beta)
                self.system.compute_grid_gradient_fock(self.grid.points, self.grid.weights, gpot_beta, operator_beta)

    def compute(self):
        if self.system.wfn.closed_shell:
            rho = self.update_rho('full')
            sigma = self.update_sigma('full')
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                self._libxc_wrapper.compute_gga_exc_unpol(rho, sigma, edens)
            energy = self.grid.integrate(edens, rho)
            self.store_energy('libxc_%s' % self._name, energy)
        else:
            rho_both = self.update_rho('both')
            sigma_all = self.update_sigma('all')
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                self._libxc_wrapper.compute_gga_exc_pol(rho_both, sigma_all, edens)
            rho = self.update_rho('full')
            energy = self.grid.integrate(edens, rho)
            self.store_energy('libxc_%s' % self._name, energy)
        return energy


class LibXCHybridGGA(LibXCGGA):
    '''Any Hybrid GGA functional from LibXC'''
    def __init__(self, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the 'hyb_gga_' prefix.
        '''
        LibXCEnergy.__init__(self, 'hyb_gga_' + name.lower())

    def get_exx_fraction(self):
        return self._libxc_wrapper.get_hyb_exx_fraction()
