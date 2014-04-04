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


__all__ = [
    'GridGroup', 'RGridGroup', 'UGridGroup',
    'GridObservable'
]


class GridGroup(Observable):
    def __init__(self, obasis, grid, grid_terms, label='grid_group'):
        self.grid_terms = grid_terms
        self.obasis = obasis
        self.grid = grid
        Observable.__init__(self, label)

    def _get_gga(self):
        '''Flag indicating that density gradients are used'''
        return any([grid_term.gga for grid_term in self.grid_terms])

    gga = property(_get_gga)

    def _get_potentials(self, cache):
        raise NotImplementedError

    def _update_rho(self, cache, select):
        rho, new = cache.load('rho_%s' % select, alloc=self.grid.size)
        if new:
            dm = cache['dm_%s' % select]
            self.obasis.compute_grid_density_dm(dm, self.grid.points, rho)
        return rho

    def _update_grad(self, cache, select):
        grad_rho, new = cache.load('grad_rho_%s' % select, alloc=(self.grid.size, 3))
        if new:
            dm = cache['dm_%s' % select]
            self.obasis.compute_grid_gradient_dm(dm, self.grid.points, grad_rho)
        return grad_rho

    def _update_grid_data(self, cache):
        raise NotImplementedError

    def compute(self, cache):
        # compute stuff on the grid that the grid_observables may use
        self._update_grid_data(cache)

        # compute energy terms and sum up
        result = 0.0
        for grid_term in self.grid_terms:
            energy = grid_term.compute(cache, self.grid)
            cache['energy_%s' % grid_term.label] = energy
            result += energy
        return result

    def add_fock(self, cache, *focks):
        # Get the potentials. If they are not yet evaluated, some computations
        # are needed.
        dpots, gpots, new = self._get_potentials(cache)

        if new:
            # compute stuff on the grid that the grid_observables may use
            self._update_grid_data(cache)

            # Collect the total potentials.
            args = dpots + gpots
            for grid_term in self.grid_terms:
                grid_term.add_pot(cache, self.grid, *args)

        for ichannel in xrange(len(focks)):
            # d = density
            self.obasis.compute_grid_density_fock(
                self.grid.points, self.grid.weights,
                dpots[ichannel], focks[ichannel])
            if self.gga:
                self.obasis.compute_grid_gradient_fock(
                    self.grid.points, self.grid.weights,
                    gpots[ichannel], focks[ichannel])


class RGridGroup(GridGroup):
    def _get_potentials(self, cache):
        dpot, new = cache.load('dpot_total_alpha', alloc=self.grid.size)
        dpots = [dpot]
        if self.gga:
            gpot, gnew = cache.load('gpot_total_alpha', alloc=(self.grid.size, 3))
            new |= gnew
            gpots = [gpot]
        else:
            gpots = []
        if new:
            dpot[:] = 0.0
            if self.gga:
                gpot[:] = 0.0
        return dpots, gpots, new

    def _update_grid_data(self, cache):
        rho_alpha = self._update_rho(cache, 'alpha')
        rho_full, new = cache.load('rho_full', alloc=self.grid.size)
        if new:
            rho_full[:] = rho_alpha
            rho_full *= 2
        if self.gga:
            grad_rho_alpha = self._update_grad(cache, 'alpha')
            sigma_alpha, new = cache.load('sigma_alpha', alloc=self.grid.size)
            if new:
                sigma_alpha[:] = (grad_rho_alpha**2).sum(axis=1)
            grad_rho_full, new = cache.load('grad_rho_full', alloc=(self.grid.size, 3))
            if new:
                grad_rho_full[:] = grad_rho_alpha
                grad_rho_full *= 2
            sigma_full, new = cache.load('sigma_full', alloc=self.grid.size)
            if new:
                sigma_full[:] = (grad_rho_full**2).sum(axis=1)


class UGridGroup(GridGroup):
    def _get_potentials(self, cache):
        dpot_alpha, newa = cache.load('dpot_total_alpha', alloc=self.grid.size)
        dpot_beta, newb = cache.load('dpot_total_beta', alloc=self.grid.size)
        dpots = [dpot_alpha, dpot_beta]
        new = newa or newb
        if self.gga:
            gpot_alpha, gnewa = cache.load('gpot_total_alpha', alloc=(self.grid.size, 3))
            gpot_beta, gnewb = cache.load('gpot_total_beta', alloc=(self.grid.size, 3))
            new |= gnewa or gnewb
            gpots = [gpot_alpha, gpot_beta]
        else:
            gpots = []
        if new:
            dpot_alpha[:] = 0.0
            dpot_beta[:] = 0.0
            if self.gga:
                gpot_alpha[:] = 0.0
                gpot_beta[:] = 0.0
        return dpots, gpots, new

    def _update_grid_data(self, cache):
        rho_alpha = self._update_rho(cache, 'alpha')
        rho_beta = self._update_rho(cache, 'beta')
        rho_full, new = cache.load('rho_full', alloc=self.grid.size)
        if new:
            rho_full[:] = rho_alpha
            rho_full += rho_beta
        rho_both, new = cache.load('rho_both', alloc=(self.grid.size, 2))
        if new:
            rho_both[:,0] = rho_alpha
            rho_both[:,1] = rho_beta

        if self.gga:
            grad_rho_alpha = self._update_grad(cache, 'alpha')
            grad_rho_beta = self._update_grad(cache, 'beta')
            sigma_alpha, new = cache.load('sigma_alpha', alloc=self.grid.size)
            if new:
                sigma_alpha[:] = (grad_rho_alpha**2).sum(axis=1)
            sigma_beta, new = cache.load('sigma_beta', alloc=self.grid.size)
            if new:
                sigma_beta[:] = (grad_rho_beta**2).sum(axis=1)
            sigma_cross, new = cache.load('sigma_cross', alloc=self.grid.size)
            if new:
                sigma_cross[:] = (grad_rho_alpha*grad_rho_beta).sum(axis=1)
            sigma_all, new = cache.load('sigma_all', alloc=(self.grid.size, 3))
            if new:
                sigma_all[:,0] = sigma_alpha
                sigma_all[:,1] = sigma_cross
                sigma_all[:,2] = sigma_beta


class GridObservable(object):
    gga = False

    def __init__(self, label):
        self.label = label

    def compute(self, cache, grid):
        raise NotImplementedError

    def add_pot(self, cache, grid, *args):
        raise NotImplementedError
