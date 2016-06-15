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
'''Interface to LDA, GGA and hybrid functionals from LibXC'''


import numpy as np

from horton.log import log, timer
from horton.utils import doc_inherit
from horton.meanfield.gridgroup import GridObservable
from horton.meanfield.cext import RLibXCWrapper, ULibXCWrapper


__all__ = [
    'LibXCEnergy',
    'RLibXCLDA', 'ULibXCLDA',
    'RLibXCGGA', 'ULibXCGGA',
    'RLibXCHybridGGA', 'ULibXCHybridGGA',
]


class LibXCEnergy(GridObservable):
    '''Base class for LibXC functionals'''

    prefix = None
    LibXCWrapper = None

    def __init__(self, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the ``lda_``,
                ``gga_`` or ``hyb_gga_`` prefix. (The type of functional is
                determined by the subclass.)
        '''
        name = '%s_%s' % (self.prefix, name)
        self._name = name
        self._libxc_wrapper = self.LibXCWrapper(name)
        log.cite('marques2012', 'using LibXC, the library of exchange and correlation functionals')
        GridObservable.__init__(self, 'libxc_%s' % name)



class RLibXCLDA(LibXCEnergy):
    '''Any LDA functional from LibXC for restricted wavefunctions'''

    prefix = 'lda'
    LibXCWrapper = RLibXCWrapper

    @timer.with_section('LDA edens')
    @doc_inherit(LibXCEnergy)
    def compute_energy(self, cache, grid):
        # LibXC expects the following input:
        #   - total density
        # LibXC computes:
        #   - the energy density per electron.
        rho_full = cache['rho_full']
        edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
        if new:
            self._libxc_wrapper.compute_lda_exc(rho_full, edens)
        return grid.integrate(edens, rho_full)

    @timer.with_section('LDA pot')
    @doc_inherit(LibXCEnergy)
    def add_pot(self, cache, grid, dpot_alpha):
        # LibXC expects the following input:
        #   - total density
        # LibXC computes:
        #   - the potential for the alpha electrons.
        pot, new = cache.load('pot_libxc_%s_alpha' % self._name, alloc=grid.size)
        if new:
            self._libxc_wrapper.compute_lda_vxc(cache['rho_full'], pot)
        dpot_alpha += pot


class ULibXCLDA(LibXCEnergy):
    '''Any LDA functional from LibXC for unrestricted wavefunctions'''
    prefix = 'lda'
    LibXCWrapper = ULibXCWrapper

    @timer.with_section('LDA edens')
    @doc_inherit(LibXCEnergy)
    def compute_energy(self, cache, grid):
        # LibXC expects the following input:
        #   - alpha density
        #   - beta density
        # LibXC computes:
        #   - the energy density per electron.

        # In case of spin-polarized computations, alpha and beta densities
        # go in and the 'total' energy density comes out.
        edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
        if new:
            self._libxc_wrapper.compute_lda_exc(cache['rho_both'], edens)
        return grid.integrate(edens, cache['rho_full'])

    @timer.with_section('LDA pot')
    @doc_inherit(LibXCEnergy)
    def add_pot(self, cache, grid, dpot_alpha, dpot_beta):
        # LibXC expects the following input:
        #   - alpha density
        #   - beta density
        # LibXC computes:
        #   - potential for the alpha electrons
        #   - potential for the beta electrons
        pot_both, new = cache.load('pot_libxc_%s_both' % self._name, alloc=(grid.size, 2))
        if new:
            self._libxc_wrapper.compute_lda_vxc(cache['rho_both'], pot_both)
        dpot_alpha += pot_both[:,0]
        dpot_beta += pot_both[:,1]


class RLibXCGGA(LibXCEnergy):
    gga = True
    prefix = 'gga'
    LibXCWrapper = RLibXCWrapper

    @timer.with_section('GGA edens')
    @doc_inherit(LibXCEnergy)
    def compute_energy(self, cache, grid):
        # LibXC expects the following input:
        #   - total density
        #   - norm squared of the gradient of the total density
        # LibXC computes:
        #   - energy density per electron
        rho_full = cache['rho_full']
        edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
        if new:
            sigma_full = cache['sigma_full']
            self._libxc_wrapper.compute_gga_exc(rho_full, sigma_full, edens)
        return grid.integrate(edens, rho_full)

    @timer.with_section('GGA pot')
    @doc_inherit(LibXCEnergy)
    def add_pot(self, cache, grid, dpot_alpha, gpot_alpha):
        # LibXC expects the following input:
        #   - total density
        #   - norm of the gradient of the total density
        # LibXC computes:
        #   - the derivative of the energy towards the alpha density.
        #   - the derivative of the energy towards the norm squared of the alpha density.
        dpot, newd = cache.load('dpot_libxc_%s_alpha' % self._name, alloc=grid.size)
        spot, news = cache.load('spot_libxc_%s_alpha' % self._name, alloc=grid.size)
        if newd or news:
            rho_full = cache['rho_full']
            sigma_full = cache['sigma_full']
            self._libxc_wrapper.compute_gga_vxc(rho_full, sigma_full, dpot, spot)

        gpot, new = cache.load('gpot_libxc_%s_alpha' % self._name, alloc=(grid.size,3))
        if new:
            grad_rho = cache['grad_rho_full']
            np.multiply(grad_rho, spot.reshape(-1,1), out=gpot)
            gpot *= 2

        dpot_alpha += dpot
        gpot_alpha += gpot


class ULibXCGGA(LibXCEnergy):
    gga = True
    prefix = 'gga'
    LibXCWrapper = ULibXCWrapper

    @timer.with_section('GGA edens')
    @doc_inherit(LibXCEnergy)
    def compute_energy(self, cache, grid):
        # LibXC expects the following input:
        #   - alpha density
        #   - beta density
        #   - norm squared of the gradient of the alpha density
        #   - dot product of the gradient of the alpha and beta densities
        #   - norm squared of the gradient of the beta density
        # LibXC computes:
        #   - energy density per electron
        edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
        if new:
            rho_both = cache['rho_both']
            sigma_all = cache['sigma_all']
            self._libxc_wrapper.compute_gga_exc(rho_both, sigma_all, edens)
        rho_full = cache['rho_full']
        return grid.integrate(edens, rho_full)

    @doc_inherit(LibXCEnergy)
    @timer.with_section('GGA pot')
    def add_pot(self, cache, grid, dpot_alpha, dpot_beta, gpot_alpha, gpot_beta):
        # LibXC expects the following input:
        #   - alpha density
        #   - beta density
        #   - norm squared of the gradient of the alpha density
        #   - dot product of the gradient of the alpha and beta densities
        #   - norm squared of the gradient of the beta density
        # LibXC computes:
        #   - the derivative of the energy towards the alpha density.
        #   - the derivative of the energy towards the beta density.
        #   - the derivative of the energy towards the norm squared of the alpha density.
        #   - the derivative of the energy towards the dot product of the alpha and beta densities.
        #   - the derivative of the energy towards the norm squared of the beta density.
        dpot_both, newd = cache.load('dpot_libxc_%s_both' % self._name, alloc=(grid.size, 2))
        spot_all, newt = cache.load('spot_libxc_%s_all' % self._name, alloc=(grid.size, 3))
        if newd or newt:
            rho_both = cache['rho_both']
            sigma_all = cache['sigma_all']
            self._libxc_wrapper.compute_gga_vxc(rho_both, sigma_all, dpot_both, spot_all)

        gpot_xc_alpha, new = cache.load('gpot_libxc_%s_alpha' % self._name, alloc=(grid.size,3))
        if new:
            gpot_xc_alpha[:] = (2*spot_all[:,0].reshape(-1,1))*cache['grad_rho_alpha']
            gpot_xc_alpha[:] += (spot_all[:,1].reshape(-1,1))*cache['grad_rho_beta']

        gpot_xc_beta, new = cache.load('gpot_libxc_%s_beta' % self._name, alloc=(grid.size,3))
        if new:
            gpot_xc_beta[:] = (2*spot_all[:,2].reshape(-1,1))*cache['grad_rho_beta']
            gpot_xc_beta[:] += (spot_all[:,1].reshape(-1,1))*cache['grad_rho_alpha']

        dpot_alpha += dpot_both[:,0]
        dpot_beta += dpot_both[:,1]
        gpot_alpha += gpot_xc_alpha
        gpot_beta += gpot_xc_beta


class RLibXCHybridGGA(RLibXCGGA):
    '''Any Hybrid GGA functional from LibXC for restricted wavefunctions'''
    prefix = 'hyb_gga'

    def get_exx_fraction(self):
        '''Return the fraction of Hartree-Fock exchange for this functional'''
        return self._libxc_wrapper.get_hyb_exx_fraction()


class ULibXCHybridGGA(ULibXCGGA):
    '''Any Hybrid GGA functional from LibXC for unrestricted wavefunctions'''
    prefix = 'hyb_gga'

    def get_exx_fraction(self):
        '''Return the fraction of Hartree-Fock exchange for this functional'''
        return self._libxc_wrapper.get_hyb_exx_fraction()
