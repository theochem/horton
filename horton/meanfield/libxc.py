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
from horton.meanfield.gridgroup import GridObservable
from horton.meanfield.cext import LibXCWrapper
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = ['LibXCLDA', 'LibXCGGA', 'LibXCHybridGGA']



class LibXCEnergy(GridObservable):
    def __init__(self, wfn, prefix, name):
        name = '%s_%s' % (prefix, name)
        self._name = name
        self._libxc_wrapper = LibXCWrapper(name)
        self._wfn = wfn
        log.cite('marques2012', 'using LibXC, the library of exchange and correlation functionals')
        GridObservable.__init__(self, 'libxc_%s' % name)


class LibXCLDA(LibXCEnergy):
    '''Any LDA functional from LibXC'''

    def __init__(self, wfn, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the ``lda_``
                prefix.
        '''
        LibXCEnergy.__init__(self, wfn, 'lda', name)

    @timer.with_section('LDA edens')
    def compute_energy(self, cache, grid):
        if isinstance(self._wfn, RestrictedWFN):
            # In the unpolarized case, libxc expects the total density as input
            # and returns the energy density per particle.
            rho_full = cache['rho_full']
            edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
            if new:
                self._libxc_wrapper.compute_lda_exc_unpol(rho_full, edens)
            return grid.integrate(edens, rho_full)
        else:
            # In case of spin-polarized computations, alpha and beta densities
            # go in and the 'total' energy density comes out.
            edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
            if new:
                rho_both = cache['rho_both']
                self._libxc_wrapper.compute_lda_exc_pol(rho_both, edens)

            rho_full = cache['rho_full']
            return grid.integrate(edens, rho_full)

    @timer.with_section('LDA pot')
    def compute_pot(self, cache, grid):
        if isinstance(self._wfn, RestrictedWFN):
            # In the closed-shell case, libxc expects the total density as input
            # and returns the potential for the alpha electrons.
            pot, new = cache.load('pot_libxc_%s_alpha' % self._name, alloc=grid.size)
            if new:
                rho_full = cache['rho_full']
                self._libxc_wrapper.compute_lda_vxc_unpol(rho_full, pot)

            cache['dpot_total_alpha'] += pot
        else:
            # In case of spin-polarized computations, alpha and beta densities
            # go in and the alpha and beta potentials come out.
            pot_both, new = cache.load('pot_libxc_%s_both' % self._name, alloc=(grid.size, 2))
            if new:
                rho_both = cache['rho_both']
                self._libxc_wrapper.compute_lda_vxc_pol(rho_both, pot_both)

            cache['dpot_total_alpha'] += pot_both[:,0]
            cache['dpot_total_beta'] += pot_both[:,1]


class LibXCGGA(LibXCEnergy):
    gga = True

    '''Any GGA functional from LibXC'''
    def __init__(self, wfn, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the ``gga_``
                prefix.
        '''
        LibXCEnergy.__init__(self, wfn, 'gga', name)

    @timer.with_section('GGA edens')
    def compute_energy(self, cache, grid):
        if isinstance(self._wfn, RestrictedWFN):
            rho_full = cache['rho_full']
            edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
            if new:
                sigma_full = cache['sigma_full']
                self._libxc_wrapper.compute_gga_exc_unpol(rho_full, sigma_full, edens)
            return grid.integrate(edens, rho_full)
        else:
            rho_both = cache['rho_both']
            edens, new = cache.load('edens_libxc_%s_full' % self._name, alloc=grid.size)
            if new:
                sigma_all = cache['sigma_all']
                self._libxc_wrapper.compute_gga_exc_pol(rho_both, sigma_all, edens)
            rho_full = cache['rho_full']
            return grid.integrate(edens, rho_full)

    @timer.with_section('GGA pot')
    def compute_pot(self, cache, grid):
        if isinstance(self._wfn, RestrictedWFN):
            dpot, newd = cache.load('dpot_libxc_%s_alpha' % self._name, alloc=grid.size)
            spot, news = cache.load('spot_libxc_%s_alpha' % self._name, alloc=grid.size)
            if newd or news:
                rho_full = cache['rho_full']
                sigma_full = cache['sigma_full']
                self._libxc_wrapper.compute_gga_vxc_unpol(rho_full, sigma_full, dpot, spot)

            gpot, new = cache.load('gpot_libxc_%s_alpha' % self._name, alloc=(grid.size,3))
            if new:
                grad_rho = cache['grad_rho_full']
                np.multiply(grad_rho, spot.reshape(-1,1), out=gpot)
                gpot *= 2

            cache['dpot_total_alpha'] += dpot
            cache['gpot_total_alpha'] += gpot
        else:
            dpot_both, newd = cache.load('dpot_libxc_%s_both' % self._name, alloc=(grid.size, 2))
            spot_all, newt = cache.load('spot_libxc_%s_all' % self._name, alloc=(grid.size, 3))
            if newd or newt:
                rho_both = cache['rho_both']
                sigma_all = cache['sigma_all']
                self._libxc_wrapper.compute_gga_vxc_pol(rho_both, sigma_all, dpot_both, spot_all)

            gpot_alpha, new = cache.load('gpot_libxc_%s_alpha' % self._name, alloc=(grid.size,3))
            if new:
                # TODO: make more efficient
                gpot_alpha[:] = (2*spot_all[:,0].reshape(-1,1))*cache['grad_rho_alpha']
                gpot_alpha[:] += (spot_all[:,1].reshape(-1,1))*cache['grad_rho_beta']

            gpot_beta, new = cache.load('gpot_libxc_%s_beta' % self._name, alloc=(grid.size,3))
            if new:
                # TODO: make more efficient
                gpot_beta[:] = (2*spot_all[:,2].reshape(-1,1))*cache['grad_rho_beta']
                gpot_beta[:] += (spot_all[:,1].reshape(-1,1))*cache['grad_rho_alpha']

            cache['dpot_total_alpha'] += dpot_both[:,0]
            cache['dpot_total_beta'] += dpot_both[:,1]
            cache['gpot_total_alpha'] += gpot_alpha
            cache['gpot_total_beta'] += gpot_beta


class LibXCHybridGGA(LibXCGGA):
    '''Any Hybrid GGA functional from LibXC'''
    def __init__(self, wfn, name):
        '''
           **Arguments:**

           name
                The name of the functional in LibXC, without the ``hyb_gga_``
                prefix.
        '''
        LibXCEnergy.__init__(self, wfn, 'hyb_gga', name)

    def get_exx_fraction(self):
        return self._libxc_wrapper.get_hyb_exx_fraction()
