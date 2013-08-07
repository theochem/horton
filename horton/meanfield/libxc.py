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
'''Interface to LDA, GGA and Hybrid functionals from LibXC'''


import numpy as np

from horton.log import log, timer
from horton.meanfield.observable import Observable
from horton.meanfield.cext import LibXCWrapper
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = ['LibXCLDA', 'LibXCGGA', 'LibXCHybridGGA']



class LibXCEnergy(Observable):
    def __init__(self, name):
        self._name = name
        self._libxc_wrapper = LibXCWrapper(name)
        log.cite('marques2012', 'using LibXC, the library of exchange and correlation functionals')
        Observable.__init__(self, 'libxc_%s' % name)

    def _update_operator(self, postpone_grid=False):
        raise NotImplementedError

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1, postpone_grid=False):
        self._update_operator(postpone_grid)
        if not postpone_grid:
            fock_alpha.iadd(self.cache.load('op_libxc_%s_alpha' % self._name), scale)
            if isinstance(self.system.wfn, UnrestrictedWFN):
                fock_beta.iadd(self.cache.load('op_libxc_%s_beta' % self._name, ), scale)


class LibXCLDA(LibXCEnergy):
    '''Any LDA functional from LibXC'''

    require_grid = True
    def __init__(self, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the ``lda_``
                prefix.
        '''
        LibXCEnergy.__init__(self, 'lda_' + name.lower())

    @timer.with_section('LDA pot')
    def _update_operator(self, postpone_grid=False):
        if isinstance(self.system.wfn, RestrictedWFN):
            # In the closed-shell case, libxc expects the total density as input
            # and returns the potential for the alpha electrons.
            pot, new = self.cache.load('pot_libxc_%s_alpha' % self._name, alloc=self.grid.size)
            if new:
                rho = self.update_rho('full')
                self._libxc_wrapper.compute_lda_vxc_unpol(rho, pot)

            self._handle_dpot(pot, postpone_grid, 'op_libxc_%s_alpha' % self._name, 'alpha')
        else:
            # In case of spin-polarized computations, alpha and beta densities
            # go in and the alpha and beta potentials come out of  it.
            pot_both, new = self.cache.load('pot_libxc_%s_both' % self._name, alloc=(self.grid.size, 2))
            if new:
                rho_both = self.update_rho('both')
                self._libxc_wrapper.compute_lda_vxc_pol(rho_both, pot_both)

            self._handle_dpot(pot_both[:,0], postpone_grid, 'op_libxc_%s_alpha' % self._name, 'alpha')
            self._handle_dpot(pot_both[:,1], postpone_grid, 'op_libxc_%s_beta' % self._name, 'beta')

    @timer.with_section('LDA edens')
    def compute(self):
        if isinstance(self.system.wfn, RestrictedWFN):
            # In the unpolarized case, libxc expects the total density as input
            # and returns the energy density per particle.
            rho = self.update_rho('full')
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                self._libxc_wrapper.compute_lda_exc_unpol(rho, edens)
            return self.grid.integrate(edens, rho)
        else:
            # In case of spin-polarized computations, alpha and beta densities
            # go in and the 'total' energy density comes out.
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                rho_both = self.update_rho('both')
                self._libxc_wrapper.compute_lda_exc_pol(rho_both, edens)

            rho = self.update_rho('full')
            return self.grid.integrate(edens, rho)


class LibXCGGA(LibXCEnergy):
    '''Any GGA functional from LibXC'''
    def __init__(self, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the ``gga_``
                prefix.
        '''
        LibXCEnergy.__init__(self, 'gga_' + name.lower())

    @timer.with_section('GGA pot')
    def _update_operator(self, postpone_grid=False):
        if isinstance(self.system.wfn, RestrictedWFN):
            dpot, newd = self.cache.load('dpot_libxc_%s_alpha' % self._name, alloc=self.grid.size)
            spot, news = self.cache.load('spot_libxc_%s_alpha' % self._name, alloc=self.grid.size)
            if newd or news:
                rho = self.update_rho('full')
                sigma = self.update_sigma('full')
                self._libxc_wrapper.compute_gga_vxc_unpol(rho, sigma, dpot, spot)

            gpot, new = self.cache.load('gpot_libxc_%s_alpha' % self._name, alloc=(self.grid.size,3))
            if new:
                grad_rho = self.update_grad_rho('full')
                np.multiply(grad_rho, spot.reshape(-1,1), out=gpot)
                gpot *= 2

            self._handle_dpot(dpot, postpone_grid, 'op_libxc_%s_alpha' % self._name, 'alpha')
            self._handle_gpot(gpot, postpone_grid, 'op_libxc_%s_alpha' % self._name, 'alpha')
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

            self._handle_dpot(dpot_both[:,0], postpone_grid, 'op_libxc_%s_alpha' % self._name, 'alpha')
            self._handle_dpot(dpot_both[:,1], postpone_grid, 'op_libxc_%s_beta' % self._name, 'beta')
            self._handle_gpot(gpot_alpha, postpone_grid, 'op_libxc_%s_alpha' % self._name, 'alpha')
            self._handle_gpot(gpot_beta, postpone_grid, 'op_libxc_%s_beta' % self._name, 'beta')

    @timer.with_section('GGA edens')
    def compute(self):
        if isinstance(self.system.wfn, RestrictedWFN):
            rho = self.update_rho('full')
            sigma = self.update_sigma('full')
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                self._libxc_wrapper.compute_gga_exc_unpol(rho, sigma, edens)
            return self.grid.integrate(edens, rho)
        else:
            rho_both = self.update_rho('both')
            sigma_all = self.update_sigma('all')
            edens, new = self.cache.load('edens_libxc_%s_full' % self._name, alloc=self.grid.size)
            if new:
                self._libxc_wrapper.compute_gga_exc_pol(rho_both, sigma_all, edens)
            rho = self.update_rho('full')
            return self.grid.integrate(edens, rho)


class LibXCHybridGGA(LibXCGGA):
    '''Any Hybrid GGA functional from LibXC'''
    def __init__(self, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the ``hyb_gga_``
                prefix.
        '''
        LibXCEnergy.__init__(self, 'hyb_gga_' + name.lower())

    def get_exx_fraction(self):
        return self._libxc_wrapper.get_hyb_exx_fraction()
