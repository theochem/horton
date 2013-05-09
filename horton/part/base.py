# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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

from horton.cache import JustOnceClass, just_once, Cache
from horton.log import log, timer
from horton.moments import get_cartesian_powers


__all__ = ['Part', 'WPart', 'CPart']


class Part(JustOnceClass):
    name = None

    def __init__(self, system, grid, local, moldens=None):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           grid
                The integration grid

           local
                Whether or not to use local (non-periodic) grids.

           **Optional arguments:**

           moldens
                The all-electron density grid data.
        '''
        JustOnceClass.__init__(self)
        self._system = system
        self._grid = grid
        self._local = local

        # Caching stuff, to avoid recomputation of earlier results
        self._cache = Cache()
        # Caching of work arrays to avoid reallocation
        if moldens is not None:
            self._cache.dump('moldens', moldens)

        # Initialize the subgrids
        if local:
            self._init_subgrids()

        # Some screen logging
        self._init_log_base()
        self._init_log_scheme()
        self._init_log_memory()

    def __getitem__(self, key):
        return self.cache.load(key)

    def _get_system(self):
        return self._system

    system = property(_get_system)

    def _get_grid(self):
        return self.get_grid()

    grid = property(_get_grid)

    def _get_local(self):
        return self._local

    local = property(_get_local)

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

    def get_grid(self, index=None):
        '''Return an integration grid

           **Optional arguments:**

           index
                The index of the atom. If not given, a grid for the entire
                system is returned. If self.local is False, a full system grid
                is always returned.
        '''
        if index is None or not self.local:
            return self._grid
        else:
            return self._subgrids[index]

    def get_moldens(self, index=None, output=None):
        self.do_moldens()
        moldens = self.cache.load('moldens')
        result = self.to_atomic_grid(index, moldens)
        if output is not None:
            output[:] = result
        return result

    def get_wcor(self, index):
        '''Return the weight corrections on a grid

           See get_grid for the meaning of the optional arguments
        '''
        raise NotImplementedError

    def _init_subgrids(self):
        raise NotImplementedError

    def _init_log_base(self):
        raise NotImplementedError

    def _init_log_scheme(self):
        raise NotImplementedError

    def _init_log_memory(self):
        if log.do_medium:
            # precompute arrays sizes for certain grids
            nbyte_global = self.grid.size*8
            nbyte_locals = np.array([self.get_grid(i).size*8 for i in xrange(self.natom)])

            # compute and report usage
            estimates = self.get_memory_estimates()
            nbyte_total = 0
            log('Coarse estimate of memory usage:')
            log('                         Label  Memory[GB]')
            log.hline()
            for label, nlocals, nglobal in estimates:
                nbyte = np.dot(nlocals, nbyte_locals) + nglobal*nbyte_global
                log('%30s  %10.3f' % (label, nbyte/1024.0**3))
                nbyte_total += nbyte
            log('%30s  %10.3f' % ('Total', nbyte_total/1024.0**3))
            log.hline()

    def get_memory_estimates(self):
        return [
            ('Atomic weights', np.ones(self.natom), 0),
            ('Promolecule', np.zeros(self.natom), 1),
            ('Working arrays', np.zeros(self.natom), 2),
        ]

    def to_atomic_grid(self, index, data):
        raise NotImplementedError

    def compute_pseudo_population(self, index):
        grid = self.get_grid(index)
        dens = self.get_moldens(index)
        at_weights = self.cache.load('at_weights', index)
        wcor = self.get_wcor(index)
        return grid.integrate(at_weights, dens, wcor)

    @just_once
    def do_moldens(self):
        raise NotImplementedError

    @just_once
    def do_partitioning(self):
        self.update_at_weights()
    do_partitioning.names = []

    def update_at_weights(self):
        '''Updates the at_weights arrays in the case (and all related arrays)'''
        raise NotImplementedError

    @just_once
    def do_populations(self):
        if log.do_medium:
            log('Computing atomic populations.')
        populations, new = self.cache.load('populations', alloc=self.system.natom)
        if new:
            self.do_partitioning()
            pseudo_populations = self.cache.load('pseudo_populations', alloc=self.system.natom)[0]
            for i in xrange(self.natom):
                pseudo_populations[i] = self.compute_pseudo_population(i)
            populations[:] = pseudo_populations
            populations += self.system.numbers - self.system.pseudo_numbers
    do_populations.names = ['populations', 'pseudo_populations']

    @just_once
    def do_charges(self):
        if log.do_medium:
            log('Computing atomic charges.')
        charges, new = self._cache.load('charges', alloc=self.system.natom)
        if new:
            self.do_populations()
            populations = self._cache.load('populations')
            charges[:] = self.system.numbers - populations
    do_charges.names = ['charges']

    @just_once
    def do_moments(self):
        if log.do_medium:
            log('Computing all sorts of AIM moments.')

        lmax = 4 # up to hexadecapoles
        cartesian_powers = get_cartesian_powers(lmax)
        cartesian_moments, new1 = self._cache.load('cartesian_moments', alloc=(self._system.natom, len(cartesian_powers)))

        radial_powers = np.arange(1, lmax+1)
        radial_moments, new2 = self._cache.load('radial_moments', alloc=(self._system.natom, len(radial_powers)))

        if new1 or new2:
            self.do_partitioning()
            for i in xrange(self._system.natom):
                # 1) Define a 'window' of the integration grid for this atom
                number = self._system.numbers[i]
                center = self._system.coordinates[i]
                grid = self.get_grid(i)

                # 2) Compute the AIM
                aim = self.get_moldens(i)*self.cache.load('at_weights', i)

                # 3) Compute weight corrections (TODO: needs to be assessed!)
                wcor = self.get_wcor(i)

                # 4) Compute Cartesian multipole moments
                counter = 0
                for nx, ny, nz in cartesian_powers:
                    if log.do_medium:
                        log('  moment %s%s%s' % ('x'*nx, 'y'*ny, 'z'*nz))
                    cartesian_moments[i, counter] = grid.integrate(aim, wcor, center=center, nx=nx, ny=ny, nz=nz, nr=0)
                    counter += 1

                # 5) Compute Radial moments
                for nr in radial_powers:
                    if log.do_medium:
                        log('  moment %s' % ('r'*nr))
                    radial_moments[i, nr-1] = grid.integrate(aim, wcor, center=center, nx=0, ny=0, nz=0, nr=nr)

    do_moments.names = ['cartesian_moments', 'radial_moments']

    def do_all(self):
        '''Computes all reasonable properties and returns a corresponding list of keys'''
        names = []
        for attr_name in dir(self):
            attr = getattr(self, attr_name)
            if callable(attr) and attr_name.startswith('do_') and attr_name != 'do_all':
                attr()
                names.extend(attr.names)
        return names


class WPart(Part):
    # TODO: add framework to evaluate AIM weights (and maybe other things) on
    # user-provided grids.

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
        Part.__init__(self, system, grid, local)


    def _init_log_base(self):
        if log.do_medium:
            log('Performing a density-based AIM analysis with a wavefunction as input.')
            log.deflist([
                ('Molecular grid', self._grid),
                ('System', self._system),
                ('Using local grids', self._local),
            ])

    def _init_subgrids(self):
        self._subgrids = self._grid.subgrids

    def get_wcor(self, index):
        return None

    @just_once
    def do_moldens(self):
        if log.do_medium:
            log('Computing densities on grids.')
        moldens, new = self.cache.load('moldens', alloc=self.grid.size)
        if new:
            self.system.compute_grid_density(self.grid.points, rhos=moldens)
    do_moldens.names = []

    def to_atomic_grid(self, index, data):
        if index is None or not self.local:
            return data
        else:
            grid = self.get_grid(index)
            return data[grid.begin:grid.end]



class CPart(Part):
    '''Base class for density partitioning schemes of cube files'''
    def __init__(self, system, grid, local, moldens, wcor_numbers, wcor_rcut_max=2.0, wcor_rcond=0.1):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           grid
                The uniform integration grid based on the cube file.

           local
                Whether or not to use local (non-periodic) grids.

           moldens
                The all-electron density grid data.

           wcor_numbers
                The list of element numbers for which weight corrections are
                needed.

           wcor_rcut_max
                The maximum cutoff sphere used for the weight corrections.

           wcor_rcond
                The regulatization strength for the weight correction equations.
        '''
        self._wcor_numbers = wcor_numbers
        self._wcor_rcut_max = wcor_rcut_max
        self._wcor_rcond = wcor_rcond
        Part.__init__(self, system, grid, local, moldens)

    def _get_wcor_numbers(self):
        return self._wcor_numbers

    wcor_numbers = property(_get_wcor_numbers)

    def _init_subgrids(self):
        # grids for non-periodic integrations
        self._subgrids = []
        for index in xrange(self.natom):
            center = self.system.coordinates[index]
            radius = self.get_cutoff_radius(index)
            self._subgrids.append(self.grid.get_window(center, radius))

    def _init_log_base(self):
        if log.do_medium:
            log('Performing a density-based AIM analysis with a cube file as input.')
            log.deflist([
                ('System', self.system),
                ('Uniform Integration Grid', self.grid),
                ('Grid shape', self.grid.shape),
                ('Using local grids', self._local),
                ('Mean spacing', '%10.5e' % (self.grid.grid_cell.volume**(1.0/3.0))),
                ('Weight corr. numbers', ' '.join(str(n) for n in self.wcor_numbers)),
                ('Weight corr. max rcut', '%10.5f' % self._wcor_rcut_max),
                ('Weight corr. rcond', '%10.5e' % self._wcor_rcond),
            ])

    def get_memory_estimates(self):
        if self.local:
            row = [('Weight corrections', np.array([n in self._wcor_numbers for n in self.system.numbers]), 0)]
        else:
            row = [('Weight corrections', np.zeros(self.natom), 1)]
        return Part.get_memory_estimates(self) + row

    def _get_wcor_low(self, label, get_funcs, index=None):
        # Get the functions
        if index is None or not self.local:
            funcs = []
            for i in xrange(self.natom):
                funcs.extend(get_funcs(i))
        else:
            funcs = get_funcs(index)

        # If no functions are collected, bork
        if len(funcs) == 0:
            return None

        grid = self.get_grid(index)

        if not self.local:
            index = None
        wcor, new = self.cache.load(label, index, alloc=grid.shape)
        if new:
            grid.compute_weight_corrections(funcs, output=wcor)
        return wcor

    def get_wcor(self, index=None):
        return self._get_wcor_low('wcor', self.get_wcor_funcs, index)

    def get_cutoff_radius(self, index):
        # The radius at which the weight function goes to zero
        raise NotImplementedError

    def get_wcor_funcs(self, index):
        raise NotImplementedError

    @just_once
    def do_moldens(self):
        moldens, new = self.cache.load('moldens', alloc=self.grid.shape)
        if new:
            raise NotImplementedError
    do_moldens.names = []

    def to_atomic_grid(self, index, data):
        if index is None or not self.local:
            return data
        else:
            grid = self.get_grid(index)
            result = grid.zeros()
            grid.extend(data, result)
            return result
