# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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

from horton.cache import JustOnceClass, just_once, Cache
from horton.log import log, timer


__all__ = ['Part', 'WPart', 'CPart']


class Part(JustOnceClass):
    def __init__(self, system, grid, moldens=None):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           grid
                The integration grid

           **Optional arguments:**

           moldens
                The all-electron density grid data.
        '''
        JustOnceClass.__init__(self)
        self._system = system
        self._grid = grid

        # Caching stuff, to avoid recomputation of earlier results
        self._cache = Cache()
        # Caching of work arrays to avoid reallocation
        self._work_cache = Cache()
        if moldens is not None:
            self._cache.dump('moldens', moldens)

        # Some screen logging
        self._init_log_base()
        self._init_log_scheme()

        # Initialize the subgrids
        self._init_subgrids()

        # If needed, prepare weight corrections for the integration on the
        # uniform grid
        with timer.section('Part wcor'):
            self._init_weight_corrections()

        # Do the essential part of the partitioning. All derived properties
        # are optional.
        with timer.section('Part weights'):
            self._init_partitioning()

        # Get rid of the work arrays
        self._work_cache.invalidate_all()

    def __getitem__(self, key):
        return self.cache.load(key)

    def _get_system(self):
        return self._system

    system = property(_get_system)

    def _get_grid(self):
        return self.get_grid()

    grid = property(_get_grid)

    def _get_natom(self):
        return self.system.natom

    natom = property(_get_natom)

    def _get_cache(self):
        return self._cache

    cache = property(_get_cache)

    def invalidate(self):
        '''Discard all cached results, e.g. because wfn changed'''
        JustOnceClass.invalidate(self)
        self.cache.invalidate_all()
        self._work_cache.invalidate_all()
        # immediately recompute the basics
        # TODO: For some schemes, the weights do not depend on the density
        # and recomputation of the atomic weights is a waste of time
        self._init_partitioning()

    def get_grid(self, index=None, periodic=True):
        '''Return an integration grid

           **Optional arguments:**

           index
                The index of the atom. If not given, a grid for the entire
                system is returned. If self.local is False, a full system grid
                is always returned.

           periodic
                Only relevant for periodic systems. By default, an integration
                grid is returned that yields integrals over one unit cell. When
                this option is set to False, the index must be provided and a
                grid is returned that covers the entire atom.
        '''
        raise NotImplementedError

    def get_moldens(self, index=None, periodic=True, output=None):
        '''Return the electron density on a grid

           See get_grid for the meaning of the optional arguments.
        '''
        raise NotImplementedError

    def get_at_weights(self, index, periodic=True, output=None):
        '''Return the atomic weight function on a grid

           See get_grid for the meaning of the optional arguments
        '''
        raise NotImplementedError

    def get_wcor(self, index, periodic=True, output=None):
        '''Return the weight corrections on a grid

           See get_grid for the meaning of the optional arguments
        '''
        raise NotImplementedError

    def _init_log_base(self):
        raise NotImplementedError

    def _init_log_scheme(self):
        raise NotImplementedError

    def _init_subgrids(self):
        raise NotImplementedError

    def _init_weight_corrections(self):
        raise NotImplementedError

    def _init_partitioning(self):
        raise NotImplementedError

    @just_once
    def do_moldens(self):
        raise NotImplementedError

    def compute_pseudo_population(self, index):
        grid = self.get_grid(index)
        dens = self.get_moldens(index)
        at_weights = self._work_cache.load('work0', grid.size, alloc=grid.shape)[0]
        self.get_at_weights(index, output=at_weights)
        wcor = self.get_wcor(index)
        return grid.integrate(at_weights, dens, wcor)

    @just_once
    def do_populations(self):
        if log.do_medium:
            log('Computing atomic populations.')
        populations, new = self.cache.load('populations', alloc=self.system.natom)
        if new:
            pseudo_populations = self.cache.load('pseudo_populations', alloc=self.system.natom)[0]
            for i in xrange(self.natom):
                pseudo_populations[i] = self.compute_pseudo_population(i)
            populations[:] = pseudo_populations
            populations += self.system.numbers - self.system.pseudo_numbers

    @just_once
    def do_charges(self):
        self.do_populations()
        if log.do_medium:
            log('Computing atomic charges.')
        charges, new = self._cache.load('charges', alloc=self.system.natom)
        if new:
            populations = self._cache.load('populations')
            charges[:] = self.system.numbers - populations

    @just_once
    def do_moments(self):
        if log.do_medium:
            log('Computing all sorts of AIM moments.')

        cartesian_powers = []
        lmax = 4 # up to hexadecapoles.
        for l in xrange(1, lmax+1):
            for nz in xrange(0, l+1):
                for ny in xrange(0, l-nz+1):
                    nx = l - ny - nz
                    cartesian_powers.append([nx, ny, nz])
        self._cache.dump('cartesian_powers', np.array(cartesian_powers))
        cartesian_moments = self._cache.load('cartesian_moments', alloc=(self._system.natom, len(cartesian_powers)))[0]

        radial_powers = np.arange(1, lmax+1)
        radial_moments = self._cache.load('radial_moments', alloc=(self._system.natom, len(radial_powers)))[0]
        self._cache.dump('radial_powers', radial_powers)

        for i in xrange(self._system.natom):
            # 1) Define a 'window' of the integration grid for this atom
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            grid = self.get_grid(i, periodic=False)

            # 2) Evaluate the non-periodic atomic weight in this window
            aim = grid.zeros()
            self.get_at_weights(i, periodic=False, output=aim)

            # 3) Extend the moldens over the window and multiply to obtain the
            #    AIM
            moldens = self.get_moldens(i, periodic=False)
            aim *= moldens
            del moldens

            # 4) Compute weight corrections (TODO: needs to be assessed!)
            wcor = self.get_wcor(i, periodic=False)

            # 5) Compute Cartesian multipoles
            counter = 0
            for nx, ny, nz in cartesian_powers:
                if log.do_medium:
                    log('  moment %s%s%s' % ('x'*nx, 'y'*ny, 'z'*nz))
                cartesian_moments[i, counter] = grid.integrate(aim, wcor, center=center, nx=nx, ny=ny, nz=nz, nr=0)
                counter += 1

            # 6) Compute Radial moments
            for nr in radial_powers:
                if log.do_medium:
                    log('  moment %s' % ('r'*nr))
                radial_moments[i, nr-1] = grid.integrate(aim, wcor, center=center, nx=0, ny=0, nz=0, nr=nr)

            del wcor
            del aim

    def do_all(self):
        '''Computes all reasonable properties and returns a corresponding list of keys'''
        self.do_populations()
        self.do_charges()
        self.do_moments()
        return ['populations', 'pseudo_populations', 'charges',
                'cartesian_powers', 'cartesian_moments', 'radial_powers',
                'radial_moments']


class WPart(Part):
    # TODO: add framework to evaluate AIM weights (and maybe other things) on
    # user-provided grids.

    name = None
    options = ['local']

    '''Base class for density partitioning schemes'''
    def __init__(self, system, grid, local=True):
        '''
           **Arguments:**

           grid
                A Molecular integration grid

           **Optional arguments:**

           local
                If ``True``: use the proper atomic grid for each AIM integral.
                If ``False``: use the entire molecular grid for each AIM integral.
                When set to ``True``, certain pairwise integrals are done with
                two atomic grids if needed.
        '''
        if local and grid.subgrids is None:
            raise ValueError('Atomic grids are discarded from molecular grid object, but are needed for local integrations.')
        self._local = local
        Part.__init__(self, system, grid)


    def _get_local(self):
        return self._local

    local = property(_get_local)

    def _init_log_base(self):
        if log.do_medium:
            log('Performing a density-based AIM analysis with a wavefunction as input.')
            log.deflist([
                ('Molecular grid', self._grid),
                ('System', self._system),
                ('Using local grids', self._local),
            ])

    def _init_subgrids(self):
        pass

    def _init_weight_corrections(self):
        pass

    def get_grid(self, index=None, periodic=True):
        # periodic gets ignored
        if index is None or not self.local:
            return self._grid
        else:
            return self._grid.subgrids[index]

    def get_moldens(self, index=None, periodic=True, output=None):
        # periodic gets ignored
        self.do_moldens()
        moldens = self.cache.load('moldens')
        grid = self.get_grid(index, periodic)
        if index is None or not self.local:
            result = moldens
        else:
            result = moldens[grid.begin:grid.end]
        if output is not None:
            output[:] = result
        return result

    def get_at_weights(self, index, periodic=True, output=None):
        # periodic gets ignored
        grid = self.get_grid(index)
        at_weights, new = self.cache.load('at_weights', index, alloc=grid.size)
        if new:
            self.compute_at_weights(index, at_weights)
        if output is not None:
            output[:] = at_weights
        return at_weights

    def get_wcor(self, index, periodic=True, output=None):
        # periodic gets ignored
        if output is not None:
            raise NotImplementedError
        return None

    def compute_at_weights(self, i0, output=None):
        raise NotImplementedError

    @just_once
    def do_moldens(self):
        if log.do_medium:
            log('Computing densities on grids.')
        moldens, new = self.cache.load('moldens', alloc=self.grid.size)
        if new:
            self.system.compute_grid_density(self.grid.points, rhos=moldens)


class CPart(Part):
    '''Base class for density partitioning schemes of cube files'''

    name = None

    def __init__(self, system, grid, moldens, store, wcor_numbers, wcor_rcut_max=2.0, wcor_rcond=0.1):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           grid
                The uniform integration grid based on the cube file.

           moldens
                The all-electron density grid data.

           store
                An instance of the class ArrayStore to store large working
                arrays.

           wcor_numbers
                The list of element numbers for which weight corrections are
                needed.

           wcor_rcut_max
                The maximum cutoff sphere used for the weight corrections.

           wcor_rcond
                The regulatization strength for the weight correction equations.
        '''
        # ArrayStore is used to avoid recomputation of huge arrays. This is not
        # always desirable due to memory constraints. Therefore the arrays
        # can be stored in a file or not stored at all. (See ArrayStore for
        # more details.) The convention is cpart is to use the store for large
        # arrays whose number scales with the system size, e.g. pro-atoms and
        # AIM densities. All other arrays are stored in the cache. This means
        # that after the initial setup of the pro-atoms, the partitioning schemes
        # must store sufficient details to recreate the proatoms when needed
        self._store = store
        self._wcor_numbers = wcor_numbers
        self._wcor_rcut_max = wcor_rcut_max
        self._wcor_rcond = wcor_rcond


        Part.__init__(self, system, grid, moldens)

    def _get_wcor_numbers(self):
        return self._wcor_numbers

    wcor_numbers = property(_get_wcor_numbers)

    def _init_log_base(self):
        if log.do_medium:
            log('Performing a density-based AIM analysis with a cube file as input.')
            log.deflist([
                ('System', self.system),
                ('Uniform Integration Grid', self.grid),
                ('Grid shape', self.grid.shape),
                ('Mean spacing', '%10.5e' % (self.grid.grid_cell.volume**(1.0/3.0))),
                ('Weight corr. numbers', ' '.join(str(n) for n in self.wcor_numbers)),
                ('Weight corr. max rcut', '%10.5f' % self._wcor_rcut_max),
                ('Weight corr. rcond', '%10.5e' % self._wcor_rcond),
            ])

    def _init_subgrids(self):
        # Windows for non-periodic integrations
        self._windows = []
        for index in xrange(self.natom):
            center = self.system.coordinates[index]
            radius = self.get_cutoff_radius(index)
            self._windows.append(self.grid.get_window(center, radius))

    def get_grid(self, index=None, periodic=True):
        if periodic:
            return self._grid
        else:
            assert index is not None
            return self._windows[index]

    def get_moldens(self, index=None, periodic=True, output=None):
        if periodic:
            result = self.cache.load('moldens')
            if output is not None:
                output[:] = result
            return result
        else:
            assert index is not None
            window = self._windows[index]
            if output is None:
                output = window.zeros()
            else:
                output[:] = 0.0
            window.extend(self.cache.load('moldens'), output)
            return output

    def get_at_weights(self, index=None, periodic=True, output=None):
        if periodic:
            # The default behavior is load the weights from the store. If this fails,
            # they must be recomputed.
            present = self._store.load(output, 'at_weights', index)
            if output is None:
                output = self.grid.zeros()
            if not present:
                self.compute_at_weights(index, output)
                self._store.dump(output, 'at_weights', index)
        else:
            assert index is not None
            window = self._windows[index]
            if output is None:
                output = window.zeros()
            self.compute_at_weights(index, output, window)
        return output

    def get_wcor(self, index=None, periodic=True, output=None):
        if periodic:
            result = self.cache.load('wcor', default=None)
        else:
            assert index is not None
            funcs = self.get_wcor_funcs(index)
            window = self._windows[index]
            result = window.compute_weight_corrections(funcs)
        if output is not None:
            output[:] = result
        return result

    def get_cutoff_radius(self, index):
        # The radius at which the weight function goes to zero
        raise NotImplementedError

    def get_wcor_funcs(self, index):
        raise NotImplementedError

    def compute_at_weights(self, index, output, window=None):
        raise NotImplementedError

    @just_once
    def do_moldens(self):
        moldens, new = self.cache.load('moldens', alloc=self.grid.size)
        if new:
            raise NotImplementedError
