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

from horton.meanfield.observable import Observable
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = ['Hartree', 'HartreeFockExchange', 'DiracExchange']


class Hartree(Observable):
    def __init__(self, label='hartree'):
        Observable.__init__(self, label)

    def _update_coulomb(self):
        '''Recompute the Coulomb operator if it has become invalid'''
        coulomb, new = self.cache.load('op_coulomb', alloc=self.system.lf.create_one_body)
        if new:
            electron_repulsion = self.system.get_electron_repulsion()
            if isinstance(self.system.wfn, RestrictedWFN):
                electron_repulsion.apply_direct(self.system.wfn.dm_alpha, coulomb)
                coulomb.iscale(2)
            else:
                electron_repulsion.apply_direct(self.system.wfn.dm_full, coulomb)

    def compute(self):
        self._update_coulomb()
        coulomb = self.cache.load('op_coulomb')
        if isinstance(self.system.wfn, RestrictedWFN):
            return coulomb.expectation_value(self.system.wfn.dm_alpha)
        else:
            return 0.5*coulomb.expectation_value(self.system.wfn.dm_full)

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        self._update_coulomb()
        coulomb = self.cache.load('op_coulomb')
        if fock_beta is None:
            # closed shell
            fock_alpha.iadd(coulomb, scale)
        else:
            # open shell
            fock_alpha.iadd(coulomb, scale)
            fock_beta.iadd(coulomb, scale)


class HartreeFockExchange(Observable):
    def __init__(self, label='exchange_hartree_fock', fraction_exchange=1.0):
        self.fraction_exchange = fraction_exchange
        Observable.__init__(self, label)

    def _update_exchange(self):
        '''Recompute the Exchange operator(s) if invalid'''
        def helper(select):
            dm = self.system.wfn.get_dm(select)
            exchange, new = self.cache.load('op_exchange_hartree_fock_%s' % select, alloc=self.system.lf.create_one_body)
            if new:
                electron_repulsion = self.system.get_electron_repulsion()
                electron_repulsion.apply_exchange(dm, exchange)

        helper('alpha')
        if isinstance(self.system.wfn, UnrestrictedWFN):
            helper('beta')

    def compute(self):
        self._update_exchange()
        if isinstance(self.system.wfn, RestrictedWFN):
            return -self.cache.load('op_exchange_hartree_fock_alpha').expectation_value(self.system.wfn.dm_alpha)
        else:
            return -0.5*self.cache.load('op_exchange_hartree_fock_alpha').expectation_value(self.system.wfn.dm_alpha) \
                   -0.5*self.cache.load('op_exchange_hartree_fock_beta').expectation_value(self.system.wfn.dm_beta)

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        self._update_exchange()
        fock_alpha.iadd(self.cache.load('op_exchange_hartree_fock_alpha'), -self.fraction_exchange*scale)
        if fock_beta is not None:
            fock_beta.iadd(self.cache.load('op_exchange_hartree_fock_beta'), -self.fraction_exchange*scale)


# TODO: Make base class for grid functionals where alpha and beta contributions are independent.
class DiracExchange(Observable):
    '''An implementation of the Dirac Exchange Functional'''

    require_grid = True
    def __init__(self, label='exchange_dirac', coeff=None):
        '''
           **Arguments:**



           **Optional arguments:**

           coeff
                The coefficient Cx in front of the Dirac exchange energy.
                It defaults to the uniform electron gas value, i.e.
                Cx = 3/4 (3/pi)^(1/3).
        '''
        if coeff is None:
            self.coeff = 3.0/4.0*(3.0/np.pi)**(1.0/3.0)
        else:
            self.coeff = coeff
        self.derived_coeff = -self.coeff*(4.0/3.0)*2**(1.0/3.0)
        Observable.__init__(self, label)

    def _update_exchange(self):
        '''Recompute the Exchange operator(s) if invalid'''
        def helper(select):
            # update grid stuff
            rho = self.update_rho(select)
            pot, new = self.cache.load('pot_exchange_dirac_%s' % select, alloc=self.grid.size)
            if new:
                pot[:] = self.derived_coeff*rho**(1.0/3.0)

            # update operator stuff
            exchange, new = self.cache.load('op_exchange_dirac_%s' % select, alloc=self.system.lf.create_one_body)
            if new:
                self.system.compute_grid_density_fock(self.grid.points, self.grid.weights, pot, exchange)

        helper('alpha')
        if isinstance(self.system.wfn, UnrestrictedWFN):
            helper('beta')

    def compute(self):
        self._update_exchange()

        def helper(select):
            pot = self.cache.load('pot_exchange_dirac_%s' % select)
            rho = self.cache.load('rho_%s' % select)
            # TODO: this integral can also be written as an expectation value
            # of the Fock operators, which is probably more efficient. However,
            # as it is now, this functional does not use density matrices at
            # all, which is also appealing.
            return self.grid.integrate(pot, rho)

        energy = helper('alpha')
        if isinstance(self.system.wfn, RestrictedWFN):
            energy *= 2
        else:
            energy += helper('beta')
        energy *= 3.0/4.0
        return energy

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        self._update_exchange()
        fock_alpha.iadd(self.cache.load('op_exchange_dirac_alpha'), scale)
        if fock_beta is not None:
            fock_beta.iadd(self.cache.load('op_exchange_dirac_beta'), scale)
