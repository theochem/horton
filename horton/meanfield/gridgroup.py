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
'''Container for observables involving numerical integration'''


from horton.meanfield.observable import Observable
from horton.utils import doc_inherit


__all__ = [
    'GridGroup', 'RGridGroup', 'UGridGroup', 'GridObservable'
]


class GridGroup(Observable):
    '''A group of terms for the effective Hamiltonian that use numerical integration'''
    def __init__(self, obasis, grid, grid_terms, label='grid_group'):
        '''
           **Arguments:**

           obasis
                The orbital basis.

           grid
                A numerical integration grid. (must have ``points`` attribute
                and ``integrate`` method.)

           grid_terms
                The contributions to the effective Hamiltonian. This must be
                a list of instances of subclasses of
                :py:class:`GridObservable`.

           **Optional arguments:**

           label
                A label for the group.
        '''
        self.grid_terms = grid_terms
        self.obasis = obasis
        self.grid = grid
        Observable.__init__(self, label)

    def _get_gga(self):
        '''Flag indicating that density gradients are used'''
        return any([grid_term.gga for grid_term in self.grid_terms])

    gga = property(_get_gga)

    def _get_potentials(self, cache):
        '''Get list of output arrays passed to ```GridObservable.add_pot```.

           **Arguments:**

           cache
                An instance of Cache, used to store intermediate results.
        '''
        raise NotImplementedError

    def _update_rho(self, cache, select):
        '''Recompute a density when not present in the cache.

           **Arguments:**

           cache
                An instance of Cache, used to store intermediate results.

           select
                'alpha' or 'beta'.
        '''
        rho, new = cache.load('rho_%s' % select, alloc=self.grid.size)
        if new:
            dm = cache['dm_%s' % select]
            self.obasis.compute_grid_density_dm(dm, self.grid.points, rho)
        return rho

    def _update_grad(self, cache, select):
        '''Recompute a density gradient when not present in the cache.

           **Arguments:**

           cache
                An instance of Cache, used to store intermediate results.

           select
                'alpha' or 'beta'.
        '''
        grad_rho, new = cache.load('grad_rho_%s' % select, alloc=(self.grid.size, 3))
        if new:
            dm = cache['dm_%s' % select]
            self.obasis.compute_grid_gradient_dm(dm, self.grid.points, grad_rho)
        return grad_rho

    def _update_grid_data(self, cache):
        '''Compute all grid data used as input for GridObservable instances

           **Arguments:**

           cache
                An instance of Cache, used to store intermediate results.
        '''
        raise NotImplementedError

    def compute_energy(self, cache):
        '''Compute the sum of the expectation values.

           **Arguments:**

           cache
                An instance of Cache, used to store intermediate results.

           This method basically dispatches the work to all ``GridObservable``
           instances in ``self.grid_terms``.
        '''
        # compute stuff on the grid that the grid_observables may use
        self._update_grid_data(cache)

        # compute energy terms and sum up
        result = 0.0
        for grid_term in self.grid_terms:
            energy = grid_term.compute_energy(cache, self.grid)
            cache['energy_%s' % grid_term.label] = energy
            result += energy
        return result

    def add_fock(self, cache, *focks):
        '''Add contributions to the Fock matrix

           **Arguments:**

           cache
                An instance of Cache, used to store intermediate results.

           This method basically dispatches the work to all ``GridObservable``
           instances in ``self.grid_terms``.
        '''
        # Get the potentials. If they are not yet evaluated, some computations
        # are needed.
        dpots, gpots, new = self._get_potentials(cache)

        if new:
            # compute stuff on the grid that the grid_observables may use
            self._update_grid_data(cache)

            # Collect the total potentials.
            for grid_term in self.grid_terms:
                if grid_term.gga:
                    grid_term.add_pot(cache, self.grid, *(dpots + gpots))
                else:
                    grid_term.add_pot(cache, self.grid, *dpots)

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
    '''GridGroup for restricted wavefunctions.

       When the ``compute`` and ``add_pot`` methods of
       :py:class:`GridObservable` instances is called, the following functions
       are pre-computed in the integration grid and stored in the cache:

       **When LDA and/or GGA functionals are used:**

       rho_alpha
            The alpha electron density.

       rho_full
            The spin-summed electron density.

       **When LDA and/or GGA functionals are used:**

       grad_rho_alpha
            The gradient of the alpha electron density.

       grad_rho_full
            The gradient of the spin-summed electron density.

       sigma_alpha
            The norm-squared of the gradient of the alpha electron density.

       sigma_full
            The norm-squared of the gradient of the spin-summed electron density.
    '''

    @doc_inherit(GridGroup)
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

    @doc_inherit(GridGroup)
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
    '''GridGroup for unrestricted wavefunctions.

       When the ``compute`` and ``add_pot`` methods of
       :py:class:`GridObservable` instances is called, the following functions
       are pre-computed in the integration grid and stored in the cache:

       **When LDA and/or GGA functionals are used:**

       rho_alpha
            The alpha electron density.

       rho_beta
            The beta electron density.

       rho_full
            The spin-summed electron density.

       rho_both
            An array with alpha and beta electron densities. Shape=(grid.size,
            2). This is mostly useful for LibXC.

       **When LDA and/or GGA functionals are used:**

       grad_rho_alpha
            The gradient of the alpha electron density.

       grad_rho_beta
            The gradient of the alpha electron density.

       sigma_alpha
            The norm-squared of the gradient of the alpha electron density.

       sigma_cross
            The dot product of the gradient of alpha and beta electron
            densities.

       sigma_beta
            The norm-squared of the gradient of the beta electron density.

       sigma_all
            An array with all three sigma quantities combined. Shape=(grid.size,
            3). This is mostly useful for LibXC
    '''
    @doc_inherit(GridGroup)
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

    @doc_inherit(GridGroup)
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
    '''Base class for contributions to the GridGroup object'''
    gga = False

    def __init__(self, label):
        '''
           **Arguments:**

           label
                A unique label for this contribution
        '''
        self.label = label

    def compute_energy(self, cache, grid):
        '''Compute the expectation value using numerical integration

           **Arguments:**

           cache
                A Cache instance used to share intermediate results between
                the ``compute`` and ``add_pot`` methods. This cache will also
                contain pre-computed functions evaluate on the grid. See
                :py:class:`RGridGroup` and :py:class:`UGridGroup` for more
                details.

           grid
                A numerical integration grid
        '''
        raise NotImplementedError

    def add_pot(self, cache, grid, *args):
        '''Add the potential to the output arguments

           **Arguments:**

           cache
                A Cache instance used to share intermediate results between
                the ``compute`` and ``add_pot`` methods. This cache will also
                contain pre-computed functions evaluate on the grid. See
                :py:class:`RGridGroup` and :py:class:`UGridGroup` for more
                details.

           grid
                A numerical integration grid

           **Possible arguments:** depending on the subclass some of these may
           not be applicable.

           dpot_alpha
                The functional derivative of the expectation value toward the
                density of the alpha electron density. Shape = (grid.size,)

           dpot_beta
                The functional derivative of the expectation value toward the
                density of the beta electron density. Shape = (grid.size,)

           gpot_alpha
                The functional derivative of the expectation value toward the
                gradient of the alpha electron density. Shape = (grid.size, 3)

           gpot_beta
                The functional derivative of the expectation value toward the
                gradient of the beta electron density. Shape = (grid.size, 3)
        '''
        raise NotImplementedError
