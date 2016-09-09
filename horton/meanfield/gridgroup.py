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
"""Container for observables involving numerical integration"""


from horton.meanfield.observable import Observable
from horton.utils import doc_inherit


__all__ = [
    'DF_LEVEL_LDA', 'DF_LEVEL_GGA', 'DF_LEVEL_MGGA',
    'GridGroup', 'RGridGroup', 'UGridGroup', 'GridObservable'
]


# Define a few `levels` of density functionals. These are used to determine
# which properties need to be computed when using certain functionals. For LDA,
# only the density is needed. For GGA the density and the gradient are needed.
# For MGGA, the density, the gradient, the laplacian and the kinetic energy
# density are needed.
DF_LEVEL_LDA = 0
DF_LEVEL_GGA = 1
DF_LEVEL_MGGA = 2


class GridGroup(Observable):
    """Group of terms for the effective Hamiltonian that use numerical integration."""

    def __init__(self, obasis, grid, grid_terms, label='grid_group', density_cutoff=1e-9):
        """Initialize a GridGroup instance.

        Parameters
        ----------
        obasis : GOBasis
            The orbital basis.

        grid : IntGrid
            A numerical integration grid. (must have ``points`` attribute and
            ``integrate`` method.)

        grid_terms : list of GridObservable instances.
            The contributions to the effective Hamiltonian.

        label : str
            A label for the group.

        density_cutoff : float
            Whenever the density on a grid point falls below this threshold, all data for
            that grid point is set to zero. This is mainly relevant for functionals that
            use derivatives of the density or the orbitals, i.e. GGA and MGGA functionals.
        """
        self.grid_terms = grid_terms
        self.obasis = obasis
        self.grid = grid
        self.density_cutoff = density_cutoff
        Observable.__init__(self, label)

    def _get_df_level(self):
        """The density-functional level of this grid group.

        Returns
        -------
        df_level : int
            This can be any of the following:

            * ``DF_LEVEL_LDA``: only LDA functionals are used.
            * ``DF_LEVEL_GGA``: GGA (and LDA) functionals are used.
            * ``DF_LEVEL_MGGA``: MGGA (and LDA and/or GGA) functionals are used.
        """
        return max([grid_term.df_level for grid_term in self.grid_terms])

    df_level = property(_get_df_level)

    def _get_potentials(self, cache):
        """Get list of output arrays passed to ```GridObservable.add_pot```.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results.
        """
        raise NotImplementedError

    def _update_grid_basics(self, cache, select):
        """Recompute a density, gradient, ... when not present in the cache.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results.
        select : str
                'alpha' or 'beta'.
        """
        # Compute the density (and optionally derivatives, etc.) on all the grid
        # points.
        if self.df_level == DF_LEVEL_LDA:
            all_basics, new = cache.load('all_%s' % select, alloc=(self.grid.size, 1))
            if new:
                dm = cache['dm_%s' % select]
                self.obasis.compute_grid_density_dm(dm, self.grid.points, all_basics[:, 0])
        elif self.df_level == DF_LEVEL_GGA:
            all_basics, new = cache.load('all_%s' % select, alloc=(self.grid.size, 4))
            if new:
                dm = cache['dm_%s' % select]
                self.obasis.compute_grid_gga_dm(dm, self.grid.points, all_basics)
        elif self.df_level == DF_LEVEL_MGGA:
            all_basics, new = cache.load('all_%s' % select, alloc=(self.grid.size, 6))
            if new:
                dm = cache['dm_%s' % select]
                self.obasis.compute_grid_mgga_dm(dm, self.grid.points, all_basics)
        else:
            raise ValueError('Internal error: non-existent DF level.')

        # Prune grid data where the density is lower than the threshold
        if self.density_cutoff > 0:
            mask = all_basics[:, 0] < self.density_cutoff
            all_basics[mask, :] = 0.0
        return all_basics

    def _update_grid_data(self, cache):
        """Compute all grid data used as input for GridObservable instances.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results.
        """
        raise NotImplementedError

    def compute_energy(self, cache):
        """Compute the sum of the expectation values.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results.

        This method basically dispatches the work to all ``GridObservable`` instances in
        ``self.grid_terms``.
        """
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
        """Add contributions to the Fock matrix.

        Parameters
        ----------
        cache : Cache
            Used to store intermediate results.

        This method basically dispatches the work to all ``GridObservable`` instances in
        ``self.grid_terms``.
        """
        # Get the potentials. If they are not yet evaluated, some computations
        # are needed.
        pots, new = self._get_potentials(cache)

        if new:
            # compute stuff on the grid that the grid_observables may use
            self._update_grid_data(cache)

            # Collect the total potentials.
            for grid_term in self.grid_terms:
                grid_term.add_pot(cache, self.grid, *pots)

        for ichannel in xrange(len(focks)):
            if self.df_level == DF_LEVEL_LDA:
                self.obasis.compute_grid_density_fock(
                    self.grid.points, self.grid.weights,
                    pots[ichannel][:, 0], focks[ichannel])
            elif self.df_level == DF_LEVEL_GGA:
                self.obasis.compute_grid_gga_fock(
                    self.grid.points, self.grid.weights,
                    pots[ichannel], focks[ichannel])
            elif self.df_level == DF_LEVEL_MGGA:
                self.obasis.compute_grid_mgga_fock(
                    self.grid.points, self.grid.weights,
                    pots[ichannel], focks[ichannel])


class RGridGroup(GridGroup):
    """GridGroup for restricted wavefunctions.

    When the ``compute`` and ``add_pot`` methods of
    :py:class:`GridObservable` instances are called, the following functions
    are pre-computed on the integration grid and stored in the cache in
    contiguous arrays (as required by LibXC):

    **When LDA, GGA or MGGA functionals are used:**

    rho_full
        The spin-summed electron density.

    **When MGGA or GGA (combined with LDA) functionals are used:**

    grad_rho_full
        The spin-summed density gradient.

    sigma_full
        The norm-squared of the gradient of the spin-summed electron density.

    **When MGGA (combined with LDA/GGA) functionals are used:**

    lapl_full
        The spin-summed density Laplacian.

    tau_full
        The spin-summed kinetic energy density

    **Combined arrays, content depends on the types of functionals being
    used:**

    all_alpha
        An array with all relevant density data:

        * column  0:   alpha density
        * columns 1-3: alpha density gradient (x, y, z)
        * column  4:   alpha density Laplacian
        * column  5:   alpha kinetic energy density
    """

    @doc_inherit(GridGroup)
    def _get_potentials(self, cache):
        if self.df_level == DF_LEVEL_LDA:
            lda_pot_alpha, new = cache.load('lda_pot_total_alpha', alloc=self.grid.size)
            if new:
                lda_pot_alpha[:] = 0.0
            return (lda_pot_alpha.reshape(-1, 1),), new
        elif self.df_level == DF_LEVEL_GGA:
            gga_pot_alpha, new = cache.load('gga_pot_total_alpha', alloc=(self.grid.size, 4))
            if new:
                gga_pot_alpha[:] = 0.0
            return (gga_pot_alpha,), new
        elif self.df_level == DF_LEVEL_MGGA:
            mgga_pot_alpha, new = cache.load('mgga_pot_total_alpha', alloc=(self.grid.size, 6))
            if new:
                mgga_pot_alpha[:] = 0.0
            return (mgga_pot_alpha,), new
        else:
            raise ValueError('Internal error: non-existent DF level.')

    @doc_inherit(GridGroup)
    def _update_grid_data(self, cache):
        all_alpha = self._update_grid_basics(cache, 'alpha')
        # Compute some derived quantities
        if self.df_level >= DF_LEVEL_LDA:
            rho_full, new = cache.load('rho_full', alloc=self.grid.size)
            if new:
                rho_full[:] = 2*all_alpha[:, 0]
        if self.df_level >= DF_LEVEL_GGA:
            grad_rho_full, new = cache.load('grad_rho_full', alloc=(self.grid.size, 3))
            if new:
                grad_rho_full[:] = all_alpha[:, 1:4]
                grad_rho_full *= 2
            sigma_full, new = cache.load('sigma_full', alloc=self.grid.size)
            if new:
                sigma_full[:] = 4*(all_alpha[:, 1:4]**2).sum(axis=1)
        if self.df_level >= DF_LEVEL_MGGA:
            lapl_full, new = cache.load('lapl_full', alloc=self.grid.size)
            if new:
                lapl_full[:] = 2*all_alpha[:, 4]
            tau_full, new = cache.load('tau_full', alloc=self.grid.size)
            if new:
                tau_full[:] = 2*all_alpha[:, 5]


class UGridGroup(GridGroup):
    """GridGroup for unrestricted wavefunctions.

    When the ``compute`` and ``add_pot`` methods of
    :py:class:`GridObservable` instances is called, the following functions
    are pre-computed on the integration grid and stored in the cache in
    contiguous arrays (as required by LibXC):

    **When LDA, GGA or MGGA functionals are used:**

    rho_full
        The spin-summed electron density.

    rho_both
        An array with alpha and beta electron densities. Shape=(grid.size,
        2). This is mostly useful for LibXC.

    **When MGGA or GGA (combined with LDA) functionals are used:**

    grad_rho_full
        The spin-summed density gradient.

    sigma_all
        An array with all three sigma quantities combined. Shape=(grid.size,
        3). This is mostly useful for LibXC

    **When MGGA (combined with LDA/GGA) functionals are used:**

    lapl_both
        The Laplacian of the alpha and the beta density. Shape=(grid.size,
        2). This is mostly useful for LibXC.

    tau_both
        The alpha and beta kinetic energy density. Shape=(grid.size, 2).
        This is mostly useful for LibXC.

    **Combined arrays, content depends on the types of functionals being
    used:**

    all_alpha, all_beta
        An array with all relevant density data:

        * column  0:   alpha/beta density
        * columns 1-3: alpha/beta density gradient (x, y, z)
        * column  4:   alpha/beta density Laplacian
        * column  5:   alpha/beta kinetic energy density
    """

    @doc_inherit(GridGroup)
    def _get_potentials(self, cache):
        if self.df_level == DF_LEVEL_LDA:
            lda_pot_alpha, newa = cache.load('lda_pot_total_alpha', alloc=self.grid.size)
            if newa:
                lda_pot_alpha[:] = 0.0
            lda_pot_beta, newb = cache.load('lda_pot_total_beta', alloc=self.grid.size)
            if newb:
                lda_pot_beta[:] = 0.0
            return (lda_pot_alpha.reshape(-1, 1), lda_pot_beta.reshape(-1, 1)), (newa or newb)
        elif self.df_level == DF_LEVEL_GGA:
            gga_pot_alpha, newa = cache.load('gga_pot_total_alpha', alloc=(self.grid.size, 4))
            if newa:
                gga_pot_alpha[:] = 0.0
            gga_pot_beta, newb = cache.load('gga_pot_total_beta', alloc=(self.grid.size, 4))
            if newb:
                gga_pot_beta[:] = 0.0
            return (gga_pot_alpha, gga_pot_beta), (newa or newb)
        elif self.df_level == DF_LEVEL_MGGA:
            mgga_pot_alpha, newa = cache.load('mgga_pot_total_alpha', alloc=(self.grid.size, 6))
            if newa:
                mgga_pot_alpha[:] = 0.0
            mgga_pot_beta, newb = cache.load('mgga_pot_total_beta', alloc=(self.grid.size, 6))
            if newb:
                mgga_pot_beta[:] = 0.0
            return (mgga_pot_alpha, mgga_pot_beta), (newa or newb)
        else:
            raise ValueError('Internal error: non-existent DF level.')

    @doc_inherit(GridGroup)
    def _update_grid_data(self, cache):
        all_alpha = self._update_grid_basics(cache, 'alpha')
        all_beta = self._update_grid_basics(cache, 'beta')
        # Compute some derived quantities
        if self.df_level >= DF_LEVEL_LDA:
            rho_full, new = cache.load('rho_full', alloc=self.grid.size)
            if new:
                rho_full[:] = all_alpha[:, 0] + all_beta[:, 0]
            rho_both, new = cache.load('rho_both', alloc=(self.grid.size, 2))
            if new:
                rho_both[:, 0] = all_alpha[:, 0]
                rho_both[:, 1] = all_beta[:, 0]
        if self.df_level >= DF_LEVEL_GGA:
            grad_rho_full, new = cache.load('grad_rho_full', alloc=(self.grid.size, 3))
            if new:
                grad_rho_full[:] = all_alpha[:, 1:4]
                grad_rho_full += all_beta[:, 1:4]
            sigma_all, new = cache.load('sigma_all', alloc=(self.grid.size, 3))
            if new:
                sigma_all[:, 0] = (all_alpha[:, 1:4]**2).sum(axis=1)
                sigma_all[:, 1] = (all_alpha[:, 1:4]*all_beta[:, 1:4]).sum(axis=1)
                sigma_all[:, 2] = (all_beta[:, 1:4]**2).sum(axis=1)
        if self.df_level >= DF_LEVEL_MGGA:
            lapl_both, new = cache.load('lapl_both', alloc=(self.grid.size, 2))
            if new:
                lapl_both[:, 0] = all_alpha[:, 4]
                lapl_both[:, 1] = all_beta[:, 4]
            tau_both, new = cache.load('tau_both', alloc=(self.grid.size, 2))
            if new:
                tau_both[:, 0] = all_alpha[:, 5]
                tau_both[:, 1] = all_beta[:, 5]


class GridObservable(object):
    """Base class for contributions to the GridGroup object."""

    df_level = None

    def __init__(self, label):
        """Initialize a GridObservable.

        Parameters
        ----------
        label : str
            A unique label for this contribution.
        """
        self.label = label

    def compute_energy(self, cache, grid):
        """Compute the expectation value using numerical integration.

        Parameters
        ----------
        cache : Cache
            Used to share intermediate results between the ``compute`` and ``add_pot``
            methods. This cache will also contain pre-computed functions evaluate on the
            grid. See :py:class:`RGridGroup` and :py:class:`UGridGroup` for more details.
        grid : IntGrid
            A numerical integration grid.
        """
        raise NotImplementedError

    def add_pot(self, cache, grid, *args):
        """Add the potential to the output arguments.

        Parameters
        ----------
        cache : Cache
            Used to share intermediate results between the ``compute`` and ``add_pot``
            methods. This cache will also contain pre-computed functions evaluate on the
            grid. See :py:class:`RGridGroup` and :py:class:`UGridGroup` for more details.
        grid : IntGrid
            A numerical integration grid.
        *args : list of [np.ndarray, shape=(npoint, npot), dtype=float]
            A list of potential arrays. (Only one array for the alpha density in case of
            restricted. Two arrays, one for alpha and one for beta electrons, in case of
            unrestricted.) Each array contains `potential` data, e.g. derivatives of a
            density functional toward:

            * column 0: the density,
            * columns 1,2,3: gradient.
            * column 4: laplacian
            * column 5: kinetic energy density

            Later columns may not be present if the functional does not need them. They
            could be prexent when other terms in the effective Hamiltonian need them.
        """
        raise NotImplementedError
