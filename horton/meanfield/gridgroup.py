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
'''Container for observables involving numerical integration'''


from horton.meanfield.observable import Observable
from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = ['GridGroup', 'GridObservable']


class GridGroup(Observable):
    def __init__(self, obasis, grid, lf, wfn, grid_terms, label='grid_group'):
        self.wfn = wfn
        self.grid_terms = grid_terms
        self.obasis = obasis
        self.grid = grid
        Observable.__init__(self, lf, label)

    def _get_gga(self):
        '''Flag indicating that density gradients are used'''
        return any([grid_term.gga for grid_term in self.grid_terms])

    gga = property(_get_gga)

    def update_grid_data(self):
        def helper_rho(select):
            rho, new = self.cache.load('rho_%s' % select, alloc=self.grid.size)
            if new:
                dm = self.wfn.get_dm(select)
                self.obasis.compute_grid_density_dm(dm, self.grid.points, rho)
            return rho

        def helper_grad(select):
            grad_rho, new = self.cache.load('grad_rho_%s' % select, alloc=(self.grid.size, 3))
            if new:
                dm = self.wfn.get_dm(select)
                self.obasis.compute_grid_gradient_dm(dm, self.grid.points, grad_rho)
            return grad_rho

        if isinstance(self.wfn, RestrictedWFN):
            rho_alpha = helper_rho('alpha')
            rho_full, new = self.cache.load('rho_full', alloc=self.grid.size)
            if new:
                rho_full[:] = rho_alpha
                rho_full *= 2
            if self.gga:
                grad_rho_alpha = helper_grad('alpha')
                sigma_alpha, new = self.cache.load('sigma_alpha', alloc=self.grid.size)
                if new:
                    sigma_alpha[:] = (grad_rho_alpha**2).sum(axis=1)
                grad_rho_full, new = self.cache.load('grad_rho_full', alloc=(self.grid.size, 3))
                if new:
                    grad_rho_full[:] = grad_rho_alpha
                    grad_rho_full *= 2
                sigma_full, new = self.cache.load('sigma_full', alloc=self.grid.size)
                if new:
                    sigma_full[:] = (grad_rho_full**2).sum(axis=1)
        else:
            rho_alpha = helper_rho('alpha')
            rho_beta = helper_rho('beta')
            rho_full, new = self.cache.load('rho_full', alloc=self.grid.size)
            if new:
                rho_full[:] = rho_alpha
                rho_full += rho_beta
            rho_both, new = self.cache.load('rho_both', alloc=(self.grid.size, 2))
            if new:
                rho_both[:,0] = rho_alpha
                rho_both[:,1] = rho_beta

            if self.gga:
                grad_rho_alpha = helper_grad('alpha')
                grad_rho_beta = helper_grad('beta')
                sigma_alpha, new = self.cache.load('sigma_alpha', alloc=self.grid.size)
                if new:
                    sigma_alpha[:] = (grad_rho_alpha**2).sum(axis=1)
                sigma_beta, new = self.cache.load('sigma_beta', alloc=self.grid.size)
                if new:
                    sigma_beta[:] = (grad_rho_beta**2).sum(axis=1)
                sigma_cross, new = self.cache.load('sigma_cross', alloc=self.grid.size)
                if new:
                    sigma_cross[:] = (grad_rho_alpha*grad_rho_beta).sum(axis=1)
                sigma_all, new = self.cache.load('sigma_all', alloc=(self.grid.size, 3))
                if new:
                    sigma_all[:,0] = sigma_alpha
                    sigma_all[:,1] = sigma_cross
                    sigma_all[:,2] = sigma_beta

    def compute(self):
        # compute stuff on the grid that the grid_observables may use
        self.update_grid_data()

        # compute energy terms and sum up
        result = 0.0
        for grid_term in self.grid_terms:
            energy = grid_term.compute_energy(self.cache, self.grid)
            self.cache['energy_%s' % grid_term.label] = energy
            result += energy
        return result

    def add_fock_matrix(self, fock_alpha, fock_beta, scale=1):
        # compute stuff on the grid that the grid_observables may use
        self.update_grid_data()

        # Reset total potentials
        self.cache.load('dpot_total_alpha', alloc=self.grid.size)[0][:] = 0.0
        if self.gga:
            self.cache.load('gpot_total_alpha', alloc=(self.grid.size, 3))[0][:] = 0.0
        if isinstance(self.wfn, UnrestrictedWFN):
            self.cache.load('dpot_total_beta', alloc=self.grid.size)[0][:] = 0.0
            if self.gga:
                self.cache.load('gpot_total_beta', alloc=(self.grid.size,3))[0][:] = 0.0

        # Collect all the total potentials and turn them into contributions
        # for the fock matrix/matrices.
        for grid_term in self.grid_terms:
            grid_term.compute_pot(self.cache, self.grid)

        # d = density
        if 'dpot_total_alpha' in self.cache:
            dpot = self.cache.load('dpot_total_alpha')
            self.obasis.compute_grid_density_fock(self.grid.points, self.grid.weights, scale*dpot, fock_alpha)
        # g = gradient
        if 'gpot_total_alpha' in self.cache:
            gpot = self.cache.load('gpot_total_alpha')
            self.obasis.compute_grid_gradient_fock(self.grid.points, self.grid.weights, scale*gpot, fock_alpha)

        if isinstance(self.wfn, UnrestrictedWFN):
            # Colect potentials for beta electrons
            # d = density
            if 'dpot_total_beta' in self.cache:
                dpot = self.cache.load('dpot_total_beta')
                self.obasis.compute_grid_density_fock(self.grid.points, self.grid.weights, scale*dpot, fock_beta)
            # g = gradient
            if 'gpot_total_beta' in self.cache:
                gpot = self.cache.load('gpot_total_beta')
                self.obasis.compute_grid_gradient_fock(self.grid.points, self.grid.weights, scale*gpot, fock_beta)


class GridObservable(object):
    gga = False

    def __init__(self, label):
        self.label = label

    def compute_energy(self, cache, grid):
        raise NotImplementedError

    def compute_pot(self, cache, grid, scale=1.0):
        raise NotImplementedError
