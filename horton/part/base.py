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
'''Base classes for partitioning algorithms'''


import numpy as np

from horton.cache import JustOnceClass, just_once, Cache
from horton.log import log
from horton.moments import get_ncart_cumul, get_npure_cumul
from horton.meanfield.wfn import RestrictedWFN


__all__ = ['Part', 'WPart', 'CPart']


class Part(JustOnceClass):
    name = None
    linear = False # whether the populations are linear in the density matrix.

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

    def __clear__(self):
        self.clear()

    def clear(self):
        '''Discard all cached results, e.g. because wfn changed'''
        JustOnceClass.clear(self)
        self.cache.clear()

    def update_grid(self, grid):
        '''Specify a new grid

           **Arguments:**

           grid
                The new grid

           When the new and old grid are the same, no action is taken. When
           a really new grid is provided, the subgrids are updated and the
           cache is cleared.
        '''
        if not (grid is self._grid):
            self._grid = grid
            if self.local:
                self._init_subgrids()
            self.clear()

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

    def get_spindens(self, index=None, output=None):
        self.do_spindens()
        spindens = self.cache.load('spindens')
        result = self.to_atomic_grid(index, spindens)
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
    def do_spindens(self):
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
        populations, new = self.cache.load('populations', alloc=self.system.natom, tags='o')
        if new:
            self.do_partitioning()
            self.do_moldens()
            pseudo_populations = self.cache.load('pseudo_populations', alloc=self.system.natom, tags='o')[0]
            if log.do_medium:
                log('Computing atomic populations.')
            for i in xrange(self.natom):
                pseudo_populations[i] = self.compute_pseudo_population(i)
            populations[:] = pseudo_populations
            populations += self.system.numbers - self.system.pseudo_numbers

    @just_once
    def do_charges(self):
        charges, new = self._cache.load('charges', alloc=self.system.natom, tags='o')
        if new:
            self.do_populations()
            populations = self._cache.load('populations')
            if log.do_medium:
                log('Computing atomic charges.')
            charges[:] = self.system.numbers - populations

    @just_once
    def do_spin_charges(self):
        spin_charges, new = self._cache.load('spin_charges', alloc=self.system.natom, tags='o')
        if new:
            if isinstance(self.system.wfn, RestrictedWFN):
                spin_charges[:] = 0.0
            else:
                try:
                    self.do_spindens()
                except NotImplementedError:
                    self.cache.clear_item('spin_charges')
                    return
                self.do_partitioning()
                if log.do_medium:
                    log('Computing atomic spin charges.')
                for index in xrange(self.system.natom):
                    grid = self.get_grid(index)
                    spindens = self.get_spindens(index)
                    at_weights = self.cache.load('at_weights', index)
                    wcor = self.get_wcor(index)
                    spin_charges[index] = grid.integrate(at_weights, spindens, wcor)

    @just_once
    def do_moments(self):
        if log.do_medium:
            log('Computing cartesian and pure AIM multipoles and radial AIM moments.')

        lmax = 4 # up to hexadecapoles

        ncart = get_ncart_cumul(lmax)
        cartesian_multipoles, new1 = self._cache.load('cartesian_multipoles', alloc=(self._system.natom, ncart), tags='o')

        npure = get_npure_cumul(lmax)
        pure_multipoles, new1 = self._cache.load('pure_multipoles', alloc=(self._system.natom, npure), tags='o')

        nrad = lmax+1
        radial_moments, new2 = self._cache.load('radial_moments', alloc=(self._system.natom, nrad), tags='o')

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
                # The minus sign is present to account for the negative electron
                # charge.
                cartesian_multipoles[i] = -grid.integrate(aim, wcor, center=center, lmax=lmax, mtype=1)
                cartesian_multipoles[i, 0] += self.system.pseudo_numbers[i]

                # 5) Compute Pure multipole moments
                # The minus sign is present to account for the negative electron
                # charge.
                pure_multipoles[i] = -grid.integrate(aim, wcor, center=center, lmax=lmax, mtype=2)
                pure_multipoles[i, 0] += self.system.pseudo_numbers[i]

                # 6) Compute Radial moments
                # For the radial moments, it is not common to put a minus sign
                # for the negative electron charge.
                radial_moments[i] = grid.integrate(aim, wcor, center=center, lmax=lmax, mtype=3)

    def do_all(self):
        '''Computes all properties and return a list of their names.'''
        for attr_name in dir(self):
            attr = getattr(self, attr_name)
            if callable(attr) and attr_name.startswith('do_') and attr_name != 'do_all':
                attr()
        return list(self.cache.iterkeys(tags='o'))


class WPart(Part):
    # TODO: add framework to evaluate AIM weights (and maybe other things) on
    # user-provided grids.

    '''Base class for density partitioning schemes'''
    def __init__(self, system, grid, local=True, epsilon=0):
        '''
           **Arguments:**

           grid
                A Molecular integration grid

           **Optional arguments:**

           local
                If ``True``: use the proper atomic grid for each AIM integral.
                If ``False``: use the entire molecular grid for each AIM integral.

           epsilon
                Allow errors on the computed electron density of this magnitude
                for the sake of efficiency.
        '''
        if local and grid.subgrids is None:
            raise ValueError('Atomic grids are discarded from molecular grid object, but are needed for local integrations.')
        self._epsilon = epsilon
        Part.__init__(self, system, grid, local)

    def _get_epsilon(self):
        return self._epsilon

    epsilon = property(_get_epsilon)

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

    def _dens_helper(self, output, select='full'):
        if self.local:
            begin = 0
            pb = log.progress(self.system.natom)
            for i in xrange(self.system.natom):
                grid = self.get_grid(i)
                end = begin + grid.size
                self.system.compute_grid_density(grid.points, rhos=output[begin:end], select=select, epsilon=self.epsilon)
                begin = end
                pb()
        else:
            self.system.compute_grid_density(self.grid.points, rhos=output)

    @just_once
    def do_moldens(self):
        moldens, new = self.cache.load('moldens', alloc=self.grid.size)
        if new:
            if log.do_medium:
                log('Computing total densitiy on grids.')
            self._dens_helper(moldens)

    @just_once
    def do_spindens(self):
        spindens, new = self.cache.load('spindens', alloc=self.grid.size)
        if new:
            if log.do_medium:
                log('Computing spin densitiy on grids.')
            self._dens_helper(spindens, select='spin')

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

    @just_once
    def do_spindens(self):
        spindens, new = self.cache.load('spindens', alloc=self.grid.shape)
        if new:
            raise NotImplementedError
    do_spindens.names = []

    def to_atomic_grid(self, index, data):
        if index is None or not self.local:
            return data
        else:
            grid = self.get_grid(index)
            result = grid.zeros()
            grid.extend(data, result)
            return result
