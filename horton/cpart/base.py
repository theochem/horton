# -*- coding: utf-8 -*-
# Horton is a mol_dens Functional Theory program.
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


__all__ = ['CPart', 'CCPart']


class CPart(JustOnceClass):
    '''Base class for density partitioning schemes of cube files'''
    def __init__(self, system, ui_grid, mol_dens):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           mol_dens
                The all-electron density grid data.
        '''
        # TODO: rename mol_dens by sys_dens.
        JustOnceClass.__init__(self)

        self._system = system
        self._ui_grid = ui_grid

        # Caching stuff, to avoid recomputation of earlier results
        self._cache = Cache()
        self._cache.dump('mol_dens', mol_dens)

        # Some screen logging
        self._init_log()

        # Do the essential part of the partitioning. All derived properties
        # are optional.
        with timer.section('CPart weights'):
            self._init_at_weights()

    def __getitem__(self, key):
        return self._cache.load(key)

    def _init_log(self):
        if log.do_medium:
            log('Performing a density-based AIM analysis of a cube file.')
            log.deflist([
                ('System', self._system),
                ('Uniform Integration Grid', self._ui_grid),
                ('Grid shape', self._ui_grid.shape),
                ('Mean spacing', '%10.5e' % (self._ui_grid.grid_cell.volume**(1.0/3.0))),
            ])

    @just_once
    def _init_at_weights(self):
        raise NotImplementedError

    def _get_system(self):
        return self._system

    system = property(_get_system)

    def _get_ui_grid(self):
        return self._ui_grid

    ui_grid = property(_get_ui_grid)

    def invalidate(self):
        '''Discard all cached results, e.g. because density data changed'''
        JustOnceClass.invalidate(self)
        self._cache.invalidate()
        # immediately recompute the basics
        # TODO: For some schemes, the weights do not depend on the density
        # and recomputation of the atomic weights is a waste of time
        with timer.section('CPart weights'):
            self._init_at_weights()

    def _compute_rel_populations(self):
        '''Compute the atomic populations of the relative density'''
        result = np.zeros(self._system.natom)
        rel_mol_dens = self._cache.load('rel_mol_dens')
        for i in xrange(self._system.natom):
            at_weights = self._cache.load('at_weights', i)
            result[i] = self._ui_grid.integrate(at_weights, rel_mol_dens)
        return result

    @just_once
    def do_populations(self):
        if log.do_medium:
            log('Computing atomic populations.')
        populations, new = self._cache.load('populations', alloc=self.system.natom)
        if new:
            ref_populations = self._cache.load('ref_populations')
            populations[:] = self._compute_rel_populations() + ref_populations

    @just_once
    def do_charges(self):
        self.do_populations()
        if log.do_medium:
            log('Computing atomic charges.')
        charges, new = self._cache.load('charges', alloc=self.system.natom)
        if new:
            populations = self._cache.load('populations')
            charges[:] = self.system.numbers - populations


class CCPart(JustOnceClass):
    '''Base class for density partitioning schemes of cube files with weight corrections'''
    def __init__(self, system, ui_grid, mol_dens):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           mol_dens
                The all-electron density grid data.
        '''
        # TODO: rename mol_dens by sys_dens.
        JustOnceClass.__init__(self)

        self._system = system
        self._ui_grid = ui_grid

        # Caching stuff, to avoid recomputation of earlier results
        self._cache = Cache()
        self._cache.dump('mol_dens', mol_dens)

        # Some screen logging
        self._init_log()

        # Do the essential part of the partitioning. All derived properties
        # are optional.
        with timer.section('CCPart wcor'):
            self._init_weight_corrections()
        with timer.section('CCPart weights'):
            self._init_at_weights()

    def __getitem__(self, key):
        return self._cache.load(key)

    def _get_system(self):
        return self._system

    system = property(_get_system)

    def _get_ui_grid(self):
        return self._ui_grid

    ui_grid = property(_get_ui_grid)

    def _init_log(self):
        if log.do_medium:
            log('Performing a density-based AIM analysis of a cube file.')
            log.deflist([
                ('System', self._system),
                ('Uniform Integration Grid', self._ui_grid),
                ('Grid shape', self._ui_grid.shape),
                ('Mean spacing', '%10.5e' % (self._ui_grid.grid_cell.volume**(1.0/3.0))),
            ])

    @just_once
    def _init_weight_corrections(self):
        raise NotImplementedError

    @just_once
    def _init_at_weights(self):
        raise NotImplementedError

    def invalidate(self):
        '''Discard all cached results, e.g. because density data changed'''
        JustOnceClass.invalidate(self)
        self._cache.invalidate()
        # immediately recompute the basics
        # TODO: For some schemes, the weights do not depend on the density
        # and recomputation of the atomic weights is a waste of time
        with timer.section('CCPart weights'):
            self._init_at_weights()

    def _compute_populations(self):
        '''Compute the atomic populations'''
        result = np.zeros(self._system.natom)
        mol_dens = self._cache.load('mol_dens')
        wcor = self._cache.load('wcor')
        for i in xrange(self._system.natom):
            at_weights = self._cache.load('at_weights', i)
            result[i] = self._ui_grid.integrate(at_weights, mol_dens, wcor)
        return result

    @just_once
    def do_populations(self):
        if log.do_medium:
            log('Computing atomic populations.')
        populations, new = self._cache.load('populations', alloc=self.system.natom)
        if new:
            populations[:] = self._compute_populations()

    @just_once
    def do_charges(self):
        self.do_populations()
        if log.do_medium:
            log('Computing atomic charges.')
        charges, new = self._cache.load('charges', alloc=self.system.natom)
        if new:
            populations = self._cache.load('populations')
            charges[:] = self.system.numbers - populations
