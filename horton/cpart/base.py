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


__all__ = ['CPart']


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

        with timer.section('CPart wcor'):
            self._init_weight_corrections()
        with timer.section('CPart weights'):
            self._init_at_weights()

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

    @just_once
    def _init_weight_corrections(self):
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

    def _integrate(self, *args, **kwargs):
        wcor = self._cache.load('wcor')
        return self._ui_grid.integrate(wcor, *args, **kwargs)

    def _compute_populations(self):
        '''Compute the atomic populations'''
        result = np.zeros(self._system.natom)
        moldens = self._cache.load('moldens')
        for i in xrange(self._system.natom):
            at_weights = self._cache.load('at_weights', i)
            result[i] = self._integrate(at_weights, moldens)

        result += self.system.numbers - self.system.pseudo_numbers

        return result

    def do_all(self):
        '''Computes all reasonable properties and returns a corresponding list of keys'''
        self.do_populations()
        self.do_charges()
        self.do_dipoles()
        self.do_volumes()
        return ['populations', 'charges', 'dipoles', 'volumes']

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

    @just_once
    def do_dipoles(self):
        if log.do_medium:
            log('Computing atomic dipoles.')
        dipoles, new = self._cache.load('dipoles', alloc=(self.system.natom, 3))
        if new:
            moldens = self._cache.load('moldens')
            for index in xrange(self._system.natom):
                at_weights = self._cache.load('at_weights', index)
                center = self._system.coordinates[index]
                dipoles[index,0] = -self._integrate(at_weights, moldens, center=center, powx=1)
                dipoles[index,1] = -self._integrate(at_weights, moldens, center=center, powy=1)
                dipoles[index,2] = -self._integrate(at_weights, moldens, center=center, powz=1)

    @just_once
    def do_volumes(self):
        self.do_populations()
        if log.do_medium:
            log('Computing atomic volumes.')
        volumes, new = self._cache.load('volumes', alloc=self.system.natom,)
        if new:
            moldens = self._cache.load('moldens')
            populations = self._cache.load('populations')
            for index in xrange(self._system.natom):
                at_weights = self._cache.load('at_weights', index)
                center = self._system.coordinates[index]
                volumes[index] = self._integrate(at_weights, moldens, center=center, powr=3)/populations[index]
