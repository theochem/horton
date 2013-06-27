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
'''Base class for energy terms and other observables of the wavefunction'''


import numpy as np

from horton.log import log


__all__ = [
    'Observable',
]


class Observable(object):
    require_grid = False

    def __init__(self, label):
        self.label = label
        self._hamiltonian = None

    def set_hamiltonian(self, hamiltonian):
        if not self._hamiltonian is None:
            raise ValueError('This term is already assigned to a Hamiltonian.')
        self._hamiltonian = hamiltonian

    # The following four properties are added for convenience:

    def _get_system(self):
        '''The system for which the hamiltonian is defined.'''
        return self._hamiltonian.system

    system = property(_get_system)

    def _get_cache(self):
        '''The cache of the hamiltonian object, cleared after a change in wfn.'''
        return self._hamiltonian.cache

    cache = property(_get_cache)

    def _get_grid(self):
        '''The numerical integration grid for this hamiltonian.'''
        return self._hamiltonian.grid

    grid = property(_get_grid)

    # Generic update routines that may be useful to various base classes
    def update_rho(self, select):
        if select == 'both':
            # This is needed for libxc
            rho, new = self.cache.load('rho_both', alloc=(self.grid.size, 2))
            if new:
                rho_alpha = self.update_rho('alpha')
                rho_beta = self.update_rho('beta')
                rho[:,0] = rho_alpha
                rho[:,1] = rho_beta
        else:
            rho, new = self.cache.load('rho_%s' % select, alloc=self.grid.size)
            if new:
                self.system.compute_grid_density(self.grid.points, rhos=rho, select=select)
        return rho

    def update_grad_rho(self, select):
        grad_rho, new = self.cache.load('grad_rho_%s' % select, alloc=(self.grid.size, 3))
        if new:
            self.system.compute_grid_gradient(self.grid.points, gradrhos=grad_rho, select=select)
        return grad_rho

    def update_sigma(self, select):
        if select == 'all':
            sigma, new = self.cache.load('sigma_all', alloc=(self.grid.size, 3))
            sigma[:,0] = self.update_sigma('alpha')
            sigma[:,1] = self.update_sigma('cross')
            sigma[:,2] = self.update_sigma('beta')
        else:
            sigma, new = self.cache.load('sigma_%s' % select, alloc=self.grid.size)
            if new:
                if select == 'cross':
                    grad_rho_alpha = self.update_grad_rho('alpha')
                    grad_rho_beta = self.update_grad_rho('beta')
                    # TODO: make more efficient (with output arguments)
                    sigma[:] = (grad_rho_alpha*grad_rho_beta).sum(axis=1)
                else:
                    grad_rho = self.update_grad_rho(select)
                    # TODO: make more efficient (with output arguments)
                    sigma[:] = (grad_rho**2).sum(axis=1)
        return sigma

    def compute(self):
        raise NotImplementedError

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        raise NotImplementedError
