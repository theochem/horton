# -*- coding: utf-8 -*-
# Horton is a moldens Functional Theory program.
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


__all__ = ['CPart', 'CPart1', 'CPart2']


class CPart(JustOnceClass):
    '''Base class for density partitioning schemes of cube files'''

    name = None

    def __init__(self, system, ui_grid, moldens, smooth):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           moldens
                The all-electron density grid data.

           smooth
                When set to True, no corrections are included to integrate
                the cusps.
        '''
        JustOnceClass.__init__(self)

        self._system = system
        self._ui_grid = ui_grid
        self._smooth = smooth

        # Caching stuff, to avoid recomputation of earlier results
        self._cache = Cache()
        self._cache.dump('moldens', moldens)

        # Some screen logging
        self._init_log()

    def __getitem__(self, key):
        return self._cache.load(key)

    def _get_system(self):
        return self._system

    system = property(_get_system)

    def _get_ui_grid(self):
        return self._ui_grid

    ui_grid = property(_get_ui_grid)

    def _get_smooth(self):
        return self._smooth

    smooth = property(_get_smooth)

    @just_once
    def _init_at_weights(self):
        raise NotImplementedError

    def _init_log(self):
        if log.do_medium:
            log('Performing a density-based AIM analysis of a cube file.')
            log.deflist([
                ('System', self._system),
                ('Uniform Integration Grid', self._ui_grid),
                ('Grid shape', self._ui_grid.shape),
                ('Mean spacing', '%10.5e' % (self._ui_grid.grid_cell.volume**(1.0/3.0))),
            ])

    def invalidate(self):
        '''Discard all cached results, e.g. because density data changed'''
        JustOnceClass.invalidate(self)
        self._cache.invalidate()
        # immediately recompute the basics
        # TODO: For some schemes, the weights do not depend on the density
        # and recomputation of the atomic weights is a waste of time
        with timer.section('CPart weights'):
            self._init_at_weights()

    def _integrate(self, *args):
        return self._ui_grid.integrate(*args)

    def _compute_populations(self):
        '''Compute the atomic populations'''
        result = np.zeros(self._system.natom)
        rel_moldens = self._cache.load('moldens')
        for i in xrange(self._system.natom):
            at_weights = self._cache.load('at_weights', i)
            result[i] = self._integrate(at_weights, rel_moldens)

        nuclear_charges = self._system.props.get('nuclear_charges')
        if nuclear_charges is not None:
            result += nuclear_charges

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


class CPart1(CPart):
    '''Base class for density partitioning schemes of cube files based on correction scheme1'''
    def __init__(self, system, ui_grid, moldens, smooth):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           moldens
                The all-electron density grid data.

           smooth
                When set to True, no corrections are included to integrate
                the cusps.
        '''
        CPart.__init__(self, system, ui_grid, moldens, smooth)

        if not self.smooth:
            # Keep a copy of the original density
            self._cache.dump('orig_moldens', moldens.copy())

        # Do the essential part of the partitioning.
        with timer.section('CPart1 weights'):
            self._init_at_weights()

    def _compute_populations(self):
        result = CPart._compute_populations(self)
        if not self.smooth:
            result += self._cache.load('ref_populations')
        return result


class CPart2(CPart):
    name = None

    '''Base class for density partitioning schemes of cube files with weight corrections'''
    def __init__(self, system, ui_grid, moldens, smooth):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           moldens
                The all-electron density grid data.

           smooth
                When set to True, no corrections are included to integrate
                the cusps.
        '''
        CPart.__init__(self, system, ui_grid, moldens, smooth)

        with timer.section('CPart2 wcor'):
            self._init_weight_corrections()
        with timer.section('CPart2 weights'):
            self._init_at_weights()

    def _integrate(self, *args):
        wcor = self._cache.load('wcor')
        return self._ui_grid.integrate(wcor, *args)

    @just_once
    def _init_weight_corrections(self):
        raise NotImplementedError
