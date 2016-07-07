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
"""Interface to LDA, GGA and hybrid functionals from LibXC"""


import numpy as np

from horton.log import log, timer
from horton.utils import doc_inherit
from horton.meanfield.gridgroup import GridObservable, DF_LEVEL_LDA, \
    DF_LEVEL_GGA, DF_LEVEL_MGGA
from horton.meanfield.cext import RLibXCWrapper, ULibXCWrapper


__all__ = [
    'LibXCEnergy',
    'RLibXCLDA', 'ULibXCLDA',
    'RLibXCGGA', 'ULibXCGGA',
    'RLibXCHybridGGA', 'ULibXCHybridGGA',
    'RLibXCMGGA', 'ULibXCMGGA',
    'RLibXCHybridMGGA', 'ULibXCHybridMGGA',
]


class LibXCEnergy(GridObservable):
    """Base class for LibXC functionals."""

    prefix = None
    LibXCWrapper = None

    def __init__(self, name):
        """Initialize A LibXC functional.

        Parameters
        ----------

        name : str
            The name of the functional in LibXC, without the ``lda_``,
            ``gga_`` or ``hyb_gga_`` prefix. (The type of functional is
            determined by the subclass.)
        """
        name = '%s_%s' % (self.prefix, name)
        self._name = name
        self._libxc_wrapper = self.LibXCWrapper(name)
        log.cite('marques2012', 'using LibXC, the library of exchange and correlation functionals')
        GridObservable.__init__(self, 'libxc_%s' % name)



class RLibXCLDA(LibXCEnergy):
    """LDA functional from LibXC for restricted wavefunctions."""

    df_level = DF_LEVEL_LDA
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
    def add_pot(self, cache, grid, lda_pot_alpha):
        # LibXC expects the following input:
        #   - total density
        # LibXC computes:
        #   - the potential for the alpha electrons.
        pot, new = cache.load('pot_libxc_%s_alpha' % self._name, alloc=grid.size)
        if new:
            self._libxc_wrapper.compute_lda_vxc(cache['rho_full'], pot)
        lda_pot_alpha += pot


class ULibXCLDA(LibXCEnergy):
    """LDA functional from LibXC for unrestricted wavefunctions."""

    df_level = DF_LEVEL_LDA
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
    def add_pot(self, cache, grid, lda_pot_alpha, lda_pot_beta):
        # LibXC expects the following input:
        #   - alpha density
        #   - beta density
        # LibXC computes:
        #   - potential for the alpha electrons
        #   - potential for the beta electrons
        pot_both, new = cache.load('pot_libxc_%s_both' % self._name, alloc=(grid.size, 2))
        if new:
            self._libxc_wrapper.compute_lda_vxc(cache['rho_both'], pot_both)
        lda_pot_alpha += pot_both[:,0]
        lda_pot_beta += pot_both[:,1]


class RLibXCGGA(LibXCEnergy):
    """GGA functional from LibXC for restricted wavefunctions."""

    df_level = DF_LEVEL_GGA
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
    def add_pot(self, cache, grid, gga_pot_alpha):
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

        # Chain rule: convert derivative toward sigma into a derivative toward
        # the gradients.
        my_gga_pot_alpha, new = cache.load('gga_pot_libxc_%s_alpha' % self._name, alloc=(grid.size,4))
        if new:
            my_gga_pot_alpha[:,0] = dpot
            grad_rho = cache['grad_rho_full']
            np.multiply(grad_rho, spot.reshape(-1,1), out=my_gga_pot_alpha[:,1:4])
            my_gga_pot_alpha[:,1:4] *= 2

        # Add to the output argument
        gga_pot_alpha += my_gga_pot_alpha


class ULibXCGGA(LibXCEnergy):
    """GGA functional from LibXC for unrestricted wavefunctions."""

    df_level = DF_LEVEL_GGA
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

    @timer.with_section('GGA pot')
    @doc_inherit(LibXCEnergy)
    def add_pot(self, cache, grid, gga_pot_alpha, gga_pot_beta):
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

        # Chain rules: convert derivatives toward sigma into a derivative toward
        # the gradients.
        grad_alpha = cache['all_alpha'][:,1:4]
        grad_beta = cache['all_beta'][:,1:4]

        my_gga_pot_alpha, new = cache.load('gga_pot_libxc_%s_alpha' % self._name, alloc=(grid.size,4))
        if new:
            my_gga_pot_alpha[:,0] = dpot_both[:,0]
            my_gga_pot_alpha[:,1:4] = (2*spot_all[:,0].reshape(-1,1))*grad_alpha
            my_gga_pot_alpha[:,1:4] += (spot_all[:,1].reshape(-1,1))*grad_beta

        my_gga_pot_beta, new = cache.load('gga_pot_libxc_%s_beta' % self._name, alloc=(grid.size,4))
        if new:
            my_gga_pot_beta[:,0] = dpot_both[:,1]
            my_gga_pot_beta[:,1:4] = (2*spot_all[:,2].reshape(-1,1))*grad_beta
            my_gga_pot_beta[:,1:4] += (spot_all[:,1].reshape(-1,1))*grad_alpha

        gga_pot_alpha += my_gga_pot_alpha
        gga_pot_beta += my_gga_pot_beta


class RLibXCHybridGGA(RLibXCGGA):
    """Hybrid GGA functional from LibXC for restricted wavefunctions."""

    prefix = 'hyb_gga'

    def get_exx_fraction(self):
        """Return the fraction of Hartree-Fock exchange for this functional."""
        return self._libxc_wrapper.get_hyb_exx_fraction()


class ULibXCHybridGGA(ULibXCGGA):
    """Hybrid GGA functional from LibXC for unrestricted wavefunctions."""

    prefix = 'hyb_gga'

    def get_exx_fraction(self):
        """Return the fraction of Hartree-Fock exchange for this functional."""
        return self._libxc_wrapper.get_hyb_exx_fraction()


class RLibXCMGGA(LibXCEnergy):
    """MGGA functional from LibXC for restricted wavefunctions."""

    df_level = DF_LEVEL_MGGA
    prefix = 'mgga'
    LibXCWrapper = RLibXCWrapper

    @timer.with_section('MGGA edens')
    @doc_inherit(LibXCEnergy)
    def compute_energy(self, cache, grid):
        # LibXC expects the following input:
        #   - total density
        #   - norm squared of the gradient of the total density
        #   - the laplacian of the density
        #   - the kinetic energy density
        # LibXC computes:
        #   - energy density per electron
        rho_full = cache['rho_full']
        edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
        if new:
            sigma_full = cache['sigma_full']
            lapl_full = cache['lapl_full']
            tau_full = cache['tau_full']
            self._libxc_wrapper.compute_mgga_exc(rho_full, sigma_full, lapl_full, tau_full, edens)
        return grid.integrate(edens, rho_full)

    @timer.with_section('MGGA pot')
    @doc_inherit(LibXCEnergy)
    def add_pot(self, cache, grid, mgga_pot_alpha):
        # LibXC expects the following input:
        #   - total density
        #   - norm of the gradient of the total density
        #   - the laplacian of the density
        #   - the kinetic energy density
        # LibXC computes:
        #   - the derivative of the energy towards the alpha density.
        #   - the derivative of the energy towards the norm squared of the alpha density.
        #   - the derivative of the energy towards the laplacian of the density
        #   - the derivative of the energy towards the kinetic energy density
        dpot, newd = cache.load('dpot_libxc_%s_alpha' % self._name, alloc=grid.size)
        spot, news = cache.load('spot_libxc_%s_alpha' % self._name, alloc=grid.size)
        lpot, news = cache.load('lpot_libxc_%s_alpha' % self._name, alloc=grid.size)
        tpot, news = cache.load('tpot_libxc_%s_alpha' % self._name, alloc=grid.size)
        if newd or news:
            rho_full = cache['rho_full']
            sigma_full = cache['sigma_full']
            lapl_full = cache['lapl_full']
            tau_full = cache['tau_full']
            self._libxc_wrapper.compute_mgga_vxc(rho_full, sigma_full, lapl_full, tau_full, dpot, spot, lpot, tpot)

        # Chain rule: convert derivative toward sigma into a derivative toward
        # the gradients.
        my_mgga_pot_alpha, new = cache.load('mgga_pot_libxc_%s_alpha' % self._name, alloc=(grid.size,6))
        if new:
            my_mgga_pot_alpha[:,0] = dpot
            grad_rho = cache['grad_rho_full']
            np.multiply(grad_rho, spot.reshape(-1,1), out=my_mgga_pot_alpha[:,1:4])
            my_mgga_pot_alpha[:,1:4] *= 2
            my_mgga_pot_alpha[:,4] = lpot
            my_mgga_pot_alpha[:,5] = tpot

        # Add to the output argument
        mgga_pot_alpha += my_mgga_pot_alpha


class ULibXCMGGA(LibXCEnergy):
    """GGA functional from LibXC for unrestricted wavefunctions."""

    df_level = DF_LEVEL_MGGA
    prefix = 'mgga'
    LibXCWrapper = ULibXCWrapper

    @timer.with_section('MGGA edens')
    @doc_inherit(LibXCEnergy)
    def compute_energy(self, cache, grid):
        # LibXC expects the following input:
        #   - alpha density
        #   - beta density
        #   - norm squared of the gradient of the alpha density
        #   - dot product of the gradient of the alpha and beta densities
        #   - norm squared of the gradient of the beta density
        #   - the laplacian of the alpha density
        #   - the laplacian of the beta density
        #   - the alpha kinetic energy density
        #   - the beta kinetic energy density
        # LibXC computes:
        #   - energy density per electron
        edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
        if new:
            rho_both = cache['rho_both']
            sigma_all = cache['sigma_all']
            lapl_both = cache['lapl_both']
            tau_both = cache['tau_both']
            self._libxc_wrapper.compute_mgga_exc(rho_both, sigma_all, lapl_both, tau_both, edens)
        rho_full = cache['rho_full']
        return grid.integrate(edens, rho_full)

    @timer.with_section('MGGA pot')
    @doc_inherit(LibXCEnergy)
    def add_pot(self, cache, grid, mgga_pot_alpha, mgga_pot_beta):
        # LibXC expects the following input:
        #   - alpha density
        #   - beta density
        #   - norm squared of the gradient of the alpha density
        #   - dot product of the gradient of the alpha and beta densities
        #   - norm squared of the gradient of the beta density
        #   - the laplacian of the alpha density
        #   - the laplacian of the beta density
        #   - the alpha kinetic energy density
        #   - the beta kinetic energy density
        # LibXC computes:
        #   - the derivative of the energy towards the alpha density.
        #   - the derivative of the energy towards the beta density.
        #   - the derivative of the energy towards the norm squared of the alpha density.
        #   - the derivative of the energy towards the dot product of the alpha and beta densities.
        #   - the derivative of the energy towards the norm squared of the beta density.
        #   - the derivative of the energy towards the laplacian of the alpha density
        #   - the derivative of the energy towards the laplacian of the beta density
        #   - the derivative of the energy towards the alpha kinetic energy density
        #   - the derivative of the energy towards the beta kinetic energy density
        dpot_both, newd = cache.load('dpot_libxc_%s_both' % self._name, alloc=(grid.size, 2))
        spot_all, newt = cache.load('spot_libxc_%s_all' % self._name, alloc=(grid.size, 3))
        lpot_both, newd = cache.load('lpot_libxc_%s_both' % self._name, alloc=(grid.size, 2))
        tpot_both, newd = cache.load('tpot_libxc_%s_both' % self._name, alloc=(grid.size, 2))
        if newd or newt:
            rho_both = cache['rho_both']
            sigma_all = cache['sigma_all']
            lapl_both = cache['lapl_both']
            tau_both = cache['tau_both']
            self._libxc_wrapper.compute_mgga_vxc(rho_both, sigma_all, lapl_both, tau_both, dpot_both, spot_all, lpot_both, tpot_both)

        # Chain rules: convert derivatives toward sigma into a derivative toward
        # the gradients.
        grad_alpha = cache['all_alpha'][:,1:4]
        grad_beta = cache['all_beta'][:,1:4]

        my_mgga_pot_alpha, new = cache.load('gga_pot_libxc_%s_alpha' % self._name, alloc=(grid.size,6))
        if new:
            my_mgga_pot_alpha[:,0] = dpot_both[:,0]
            my_mgga_pot_alpha[:,1:4] = (2*spot_all[:,0].reshape(-1,1))*grad_alpha
            my_mgga_pot_alpha[:,1:4] += (spot_all[:,1].reshape(-1,1))*grad_beta
            my_mgga_pot_alpha[:,4] = lpot_both[:,0]
            my_mgga_pot_alpha[:,5] = tpot_both[:,0]

        my_mgga_pot_beta, new = cache.load('gga_pot_libxc_%s_beta' % self._name, alloc=(grid.size,6))
        if new:
            my_mgga_pot_beta[:,0] = dpot_both[:,1]
            my_mgga_pot_beta[:,1:4] = (2*spot_all[:,2].reshape(-1,1))*grad_beta
            my_mgga_pot_beta[:,1:4] += (spot_all[:,1].reshape(-1,1))*grad_alpha
            my_mgga_pot_beta[:,4] = lpot_both[:,1]
            my_mgga_pot_beta[:,5] = tpot_both[:,1]

        mgga_pot_alpha += my_mgga_pot_alpha
        mgga_pot_beta += my_mgga_pot_beta


class RLibXCHybridMGGA(RLibXCMGGA):
    """Hybrid MGGA functional from LibXC for restricted wavefunctions."""

    prefix = 'hyb_mgga'

    def get_exx_fraction(self):
        """Return the fraction of Hartree-Fock exchange for this functional."""
        return self._libxc_wrapper.get_hyb_exx_fraction()


class ULibXCHybridMGGA(ULibXCMGGA):
    """Hybrid GGA functional from LibXC for unrestricted wavefunctions."""

    prefix = 'hyb_mgga'

    def get_exx_fraction(self):
        """Return the fraction of Hartree-Fock exchange for this functional."""
        return self._libxc_wrapper.get_hyb_exx_fraction()
