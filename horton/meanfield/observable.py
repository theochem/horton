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


__all__ = [
    'Observable',
]


class Observable(object):
    require_grid = False
    exchange = False # Set to True for exhange functionals. Is needed for idiot proof option

    def __init__(self, lf, label):
        self.label = label
        self._hamiltonian = None
        self._lf = lf
        self._operator = None

    def set_hamiltonian(self, hamiltonian):
        if not self._hamiltonian is None:
            raise ValueError('This term is already assigned to a Hamiltonian.')
        self._hamiltonian = hamiltonian

    # The following four properties are added for convenience:

    def _get_system(self):
        '''The system for which the hamiltonian is defined.'''
        return self._hamiltonian.system

    system = property(_get_system)

    def _get_lf(self):
        return self._lf

    lf = property(_get_lf)

    def _get_cache(self):
        '''The cache of the hamiltonian object, cleared after a change in wfn.'''
        return self._hamiltonian.cache #FIXME: remove the hamiltonian dependency

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

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1, postpone_grid=False):
        '''Add contributions to alpha (and beta) Fock matrix(es).

           **Arguments:**

           fock_alpha
                A One-Body operator output argument for the alpha fock matrix.

           fock_alpha
                A One-Body operator output argument for the beta fock matrix.

           **Optional arguments:**

           scale
                A scale factor for this contribution

           postpone_grid
                When set to True, the grid-based observables will add their
                potential on the grid to a total for the whole hamiltonian,
                which can be converted to a fock matrix later. Note that the
                scale argument will be ignored when this option is activated.

           In the case of a closed-shell computation, the argument fock_beta is
           ``None``.
        '''
        raise NotImplementedError

    def _handle_dpot(self, dpot, postpone_grid, op_name, spin):
        '''Take care of a density potential, either make a fock contribution or collect grid data

           **Arguments:**

           dpot
                The potential grid data.

           postpone_grid
                Flag that determines how the potential is handled. When True,
                the potential grid data is added to the total for the entire
                hamiltonian. When False, the potential is converted into
                a Fock matrix and stored under the given op_name in the cache.
                When set to None, this method does nothing.

           op_name
                The name of the operator.

           spin
                'alpha', 'beta' or 'both'. This is only used if postpone_grid==True.
        '''
        if postpone_grid is True:
            # Make sure this addition is done only once
            tag = '%s_postponed_dpot' % op_name
            if tag not in self.cache:
                if spin == 'both':
                    dpot_total_alpha = self.cache.load('dpot_total_alpha', alloc=self.grid.size)[0]
                    dpot_total_alpha += dpot
                    dpot_total_beta = self.cache.load('dpot_total_beta', alloc=self.grid.size)[0]
                    dpot_total_beta += dpot
                elif spin == 'alpha' or spin == 'beta':
                    dpot_total = self.cache.load('dpot_total_%s' % spin, alloc=self.grid.size)[0]
                    dpot_total += dpot
                else:
                    raise NotImplementedError
                self.cache.dump(tag, True)
        elif postpone_grid is False:
            operator, new = self.cache.load(op_name, alloc=self.system.lf.create_one_body)
            if new:
                self.system.compute_grid_density_fock(self.grid.points, self.grid.weights, dpot, operator)
        elif postpone_grid is not None:
            raise ValueError('postpone_grid must be True, False or None')

    def _handle_gpot(self, gpot, postpone_grid, op_name, spin):
        '''Take care of a gradient potential, either make a fock contribution or collect grid data

           **Arguments:**

           gpot
                The potential grid data.

           postpone_grid
                Flag that determines how the potential is handled. When True,
                the potential grid data is added to the total for the entire
                hamiltonian. When False, the potential is converted into
                a Fock matrix and stored under the given op_name in the cache.
                When set to None, this method does nothing.

           op_name
                The name of the operator.

           spin
                'alpha', 'beta' or 'both'. This is only used if postpone_grid==True.
        '''
        if postpone_grid is True:
            # Make sure this addition is done only once
            tag = '%s_postponed_gpot' % op_name
            if tag not in self.cache:
                if spin == 'both':
                    gpot_total_alpha = self.cache.load('gpot_total_alpha', alloc=(self.grid.size,3))[0]
                    gpot_total_alpha += gpot
                    gpot_total_beta = self.cache.load('gpot_total_beta', alloc=(self.grid.size,3))[0]
                    gpot_total_beta += gpot
                elif spin == 'alpha' or spin == 'beta':
                    gpot_total = self.cache.load('gpot_total_%s' % spin, alloc=(self.grid.size,3))[0]
                    gpot_total += gpot
                else:
                    raise NotImplementedError
                self.cache.dump(tag, True)
        elif postpone_grid is False:
            operator, new = self.cache.load(op_name, alloc=self.system.lf.create_one_body)
            if new:
                self.system.compute_grid_gradient_fock(self.grid.points, self.grid.weights, gpot, operator)
        elif postpone_grid is not None:
            raise ValueError('postpone_grid must be True, False or None')
