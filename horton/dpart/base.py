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


from horton.cache import JustOnceClass, just_once, Cache
from horton.log import log, timer


__all__ = ['DPart']


class DPart(JustOnceClass):
    '''Base class for density partitioning schemes'''
    def __init__(self, molgrid, local=True):
        '''
           **Arguments:**

           molgrid
                A Molecular integration grid

           **Optional arguments:**

           local
                If ``True``: use the proper atomic grid for each AIM integral.
                If ``False``: use the entire molecular grid for each AIM integral.
                When set to ``True``, certain pairwise integrals are done with
                two atomic grids if needed.
        '''
        if local and molgrid.subgrids is None:
            raise ValueError('Atomic grids are discarded from molecular grid object, but are needed for local integrations.')

        JustOnceClass.__init__(self)

        self._molgrid = molgrid
        self._system = molgrid.system
        self._local = local

        # Caching stuff, to avoid recomputation of earlier results
        self.cache = Cache()

        # Some screen logging
        self._init_log()

        # Do the essential part of the partitioning. All derived properties
        # are optional.
        with timer.section('DPart weights'):
            self._init_at_weights()

    def __getitem__(self, key):
        return self.cache.load(key)

    def _init_log(self):
        if log.do_medium:
            log('Performing a density-based AIM analysis.')
            log.deflist([
                ('Molecular grid', self._molgrid),
                ('System', self._system),
                ('Using local grids', self._local),
            ])

    @just_once
    def _init_at_weights(self):
        raise NotImplementedError

    def _get_molgrid(self):
        return self._molgrid

    molgrid = property(_get_molgrid)

    def _get_system(self):
        return self._system

    system = property(_get_system)

    def _get_local(self):
        return self._local

    local = property(_get_local)

    def invalidate(self):
        '''Discard all cached results, e.g. because wfn changed'''
        JustOnceClass.invalidate(self)
        self.cache.invalidate()
        # immediately recompute the basics
        self._init_at_weights()

    def iter_grids(self):
        '''Iterate over the atomic grids

           **Yields:** (index, grid) pairs

           The grid may also me the same molecular grid at each iteration. This
           allows most routines to be implemented without being aware of the
           local flag. Some routines may still use the local flag to improve
           the efficiency, e.g. see do_mol_dens
        '''
        for i in xrange(self.system.natom):
            if self._local:
                yield i, self.molgrid.subgrids[i]
            else:
                yield i, self.molgrid

    @just_once
    def do_mol_dens(self):
        if log.do_medium: log('Computing densities on grids.')
        for i, grid in self.iter_grids():
            if i == 0 or self.local:
                mol_dens, new = self.cache.load('mol_dens', i, alloc=grid.size)
                if new:
                    self.system.compute_density_grid(grid.points, rhos=mol_dens)
            else:
                self.cache.dump('mol_dens', i, mol_dens)

    @just_once
    def do_populations(self):
        self.do_mol_dens()
        if log.do_medium: log('Computing atomic populations.')
        populations, new = self.cache.load('populations', alloc=self.system.natom)
        if new:
            for i, grid in self.iter_grids():
                at_weights = self.cache.load('at_weights', i)
                dens = self.cache.load('mol_dens', i)
                populations[i] = grid.integrate(at_weights, dens)

    @just_once
    def do_charges(self):
        self.do_populations()
        if log.do_medium: log('Computing atomic charges.')
        charges, new = self.cache.load('charges', alloc=self.system.natom)
        if new:
            populations = self.cache.load('populations')
            charges[:] = self.system.numbers - populations
