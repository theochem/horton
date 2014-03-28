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

from horton.meanfield.observable import Observable
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = ['Hartree', 'HartreeFockExchange', 'DiracExchange']


class Hartree(Observable):
    def __init__(self, lf, wfn, eri, label='hartree'):
        self._wfn = wfn
        self._eri = eri
        Observable.__init__(self, lf, label)

    def _update_hartree(self):
        '''Recompute the Hartree operator if it has become invalid'''
        hartree, new = self.cache.load('op_hartree',
                                        alloc=self.lf.create_one_body)
        if new:
            if isinstance(self._wfn, RestrictedWFN):
                self._eri.apply_direct(self._wfn.dm_alpha, hartree)
                hartree.iscale(2)
            else:
                self._eri.apply_direct(self._wfn.dm_full, hartree)

    def compute(self):
        self._update_hartree()
        hartree = self.cache.load('op_hartree')
        if isinstance(self._wfn, RestrictedWFN):
            return hartree.expectation_value(self._wfn.dm_alpha)
        else:
            return 0.5 * hartree.expectation_value(self._wfn.dm_full)

    def add_fock_matrix(self, fock_alpha, fock_beta,
                        scale=1, postpone_grid=False):
        self._update_hartree()
        hartree = self.cache.load('op_hartree')
        fock_alpha.iadd(hartree, scale)
        if isinstance(self._wfn, UnrestrictedWFN):
            fock_beta.iadd(hartree, scale)


class HartreeFockExchange(Observable):
    exchange = True

    def __init__(self, lf, wfn, eri,
                 label='exchange_hartree_fock', fraction_exchange=1.0):
        self._wfn = wfn
        self._eri = eri
        self.fraction_exchange = fraction_exchange
        Observable.__init__(self, lf, label)

    def _update_exchange(self):
        '''Recompute the Exchange operator(s) if invalid'''
        def helper(select):
            dm = self._wfn.get_dm(select)
            exchange, new = self.cache.load('op_exchange_hartree_fock_%s'
                                             % select, alloc=self._lf.create_one_body)
            if new:
                self._eri.apply_exchange(dm, exchange)

        helper('alpha')
        if isinstance(self._wfn, UnrestrictedWFN):
            helper('beta')

    def compute(self):
        self._update_exchange()
        xhf_fock_alpha = self.cache.load('op_exchange_hartree_fock_alpha')
        if isinstance(self._wfn, RestrictedWFN):
            return -self.fraction_exchange * xhf_fock_alpha.expectation_value(self._wfn.dm_alpha)
        else:
            xhf_fock_beta = self.cache.load('op_exchange_hartree_fock_beta')
            return -0.5 * self.fraction_exchange * xhf_fock_alpha.expectation_value(self._wfn.dm_alpha) \
                   -0.5 * self.fraction_exchange * xhf_fock_beta.expectation_value(self._wfn.dm_beta)

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1,
                         postpone_grid=False):
        self._update_exchange()
        fock_alpha.iadd(self.cache.load('op_exchange_hartree_fock_alpha'),
                         -self.fraction_exchange * scale)
        if fock_beta is not None:
            fock_beta.iadd(self.cache.load('op_exchange_hartree_fock_beta'),
                            -self.fraction_exchange * scale)


# TODO: Make base class for grid functionals where alpha and beta contributions are independent.
class DiracExchange(Observable):
    '''An implementation of the Dirac Exchange Functional'''
    require_grid = True
    exchange = True

    def __init__(self, lf, wfn, eri,
                 label='exchange_dirac', coeff=None):
        '''
           **Arguments:**



           **Optional arguments:**

           coeff
                The coefficient Cx in front of the Dirac exchange energy.
                It defaults to the uniform electron gas value, i.e.
                Cx = 3/4 (3/pi)^(1/3).
        '''
        self._wfn = wfn
        self._eri = eri
        if coeff is None:
            self.coeff = 3.0 / 4.0 * (3.0 / np.pi) ** (1.0 / 3.0)
        else:
            self.coeff = coeff
        self.derived_coeff = -self.coeff * (4.0 / 3.0) * 2 ** (1.0 / 3.0)
        Observable.__init__(self, lf, label)

    def _update_exchange(self, postpone_grid=False):
        '''Recompute the Exchange operator(s) if invalid'''
        def helper(select):
            # update grid stuff
            rho = self.update_rho(select)
            pot, new = self.cache.load('pot_exchange_dirac_%s'
                                       % select, alloc=self.grid.size)
            if new:
                pot[:] = self.derived_coeff * rho ** (1.0 / 3.0)

            # update operator stuff
            self._handle_dpot(pot, postpone_grid, 'op_exchange_dirac_%s'
                              % select, select)

        helper('alpha')
        if isinstance(self._wfn, UnrestrictedWFN):
            helper('beta')

    def compute(self):
        self._update_exchange(postpone_grid=None)

        def helper(select):
            pot = self.cache.load('pot_exchange_dirac_%s' % select)
            rho = self.cache.load('rho_%s' % select)
            return self.grid.integrate(pot, rho)

        energy = helper('alpha')
        if isinstance(self._wfn, RestrictedWFN):
            energy *= 2
        else:
            energy += helper('beta')
        energy *= 3.0 / 4.0
        return energy

    def add_fock_matrix(self, fock_alpha, fock_beta,
                        scale=1, postpone_grid=False):
        self._update_exchange(postpone_grid)
        if not postpone_grid:
            fock_alpha.iadd(self.cache.load('op_exchange_dirac_alpha'), scale)
            if fock_beta is not None:
                fock_beta.iadd(self.cache.load('op_exchange_dirac_beta'), scale)
