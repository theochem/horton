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

from horton.cpart.base import CPart
from horton.cache import just_once, Cache
from horton.log import log


__all__ = ['HirshfeldCPart', 'HirshfeldICPart']


class HirshfeldCPart(CPart):
    def __init__(self, system, ui_grid, mol_dens, proatomdb):
        '''
           **Arguments:**

           system
                The system to be partitioned.

           ui_grid
                The uniform integration grid based on the cube file.

           mol_dens
                The all-electron molecular density grid data.

           proatomdb
                The database of proatomic densities.
        '''
        self._proatomdb = proatomdb
        CPart.__init__(self, system, ui_grid, mol_dens)

    def _get_proatomdb(self):
        return self._proatomdb

    proatomdb = property(_get_proatomdb)

    def _compute_pro_atom(self, i):
        if log.do_medium:
            log('Computing pro-atom %i (%i)' % (i, self._system.natom))
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        pro_atom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
        # Construct the pro-atom
        assert new
        number = self._system.numbers[i]
        center = self._system.coordinates[i]
        spline = self._proatomdb.get_hirshfeld_proatom_fn(number)
        self._ui_grid.eval_spline(spline, center, pro_atom)
        pro_atom += 1e-100 # avoid division by zero

    def _compute_pro_molecule(self):
        # The pro-molecule is (temporarily) stored in rel_mol_dens for efficiency.
        pro_molecule, new = self._cache.load('rel_mol_dens', alloc=self._ui_grid.shape)
        if not new:
            pro_molecule[:] = 0
        for i in xrange(self._system.natom):
            pro_atom = self._cache.load('at_weights', i)
            pro_molecule += pro_atom

    def _compute_at_weights(self):
        pro_molecule = self._cache.load('rel_mol_dens')
        for i in xrange(self._system.natom):
            pro_atom = self._cache.load('at_weights', i)
            pro_atom /= pro_molecule

    def _compute_rel_dens(self):
        rel_mol_dens = self._cache.load('rel_mol_dens')
        rel_mol_dens *= -1
        rel_mol_dens += self._cache.load('mol_dens')

    @just_once
    def _init_at_weights(self):
        # A) Compute all the pro-atomic densities.
        for i in xrange(self._system.natom):
            self._compute_pro_atom(i)
        self._compute_pro_molecule()
        self._compute_at_weights()
        self._compute_rel_dens()
        # Store reference populations
        self._cache.dump('ref_populations', self._system.numbers.astype(float))


class HirshfeldICPart(HirshfeldCPart):
    def _get_isolated_atom(self, i, pop):
        # A local cache is used that only exists during _init_at_weights:
        isolated_atom, new = self._local_cache.load('isolated_atom', i, pop, alloc=self._ui_grid.shape)
        if new:
            number = self._system.numbers[i]
            center = self._system.coordinates[i]
            spline = self._proatomdb.get_hirshfeld_i_proatom_fn(number, pop)
            if log.do_medium:
                log('Computing isolated atom %i (n=%i, pop=%i)' % (i, number, pop))
            self._ui_grid.eval_spline(spline, center, isolated_atom)
        return isolated_atom

    def _compute_pro_atom(self, i, pop):
        # Pro-atoms are (temporarily) stored in at_weights for efficiency.
        pro_atom, new = self._cache.load('at_weights', i, alloc=self._ui_grid.shape)
        # Construct the pro-atom
        ipop = int(np.floor(pop))
        pro_atom[:] = self._get_isolated_atom(i, ipop)
        if float(ipop) != pop:
            x = pop - ipop
            pro_atom *= 1-x
            tmp = self._local_cache.load('tmp', alloc=self._ui_grid.shape)[0]
            tmp[:] = self._get_isolated_atom(i, ipop+1)
            tmp *= x
            pro_atom += tmp
        pro_atom += 1e-100 # avoid division by zero
        return pro_atom

    @just_once
    def _init_at_weights(self):
        self._local_cache = Cache()

        ref_populations = self._system.numbers.astype(float)
        if log.medium:
            log.hline()
            log('Iteration       Change')
            log.hline()
        counter = 0
        while True:
            # Construct pro-atoms
            for i in xrange(self._system.natom):
                self._compute_pro_atom(i, ref_populations[i])
            self._compute_pro_molecule()
            self._compute_at_weights()
            self._compute_rel_dens()
            new_populations = self._compute_rel_populations() + ref_populations
            change = abs(new_populations - ref_populations).max()
            if log.medium:
                log('%9i   %10.5e' % (counter, change))
            if change < 1e-4:
                break

            ref_populations = new_populations
            counter += 1

        if log.medium:
            log.hline()

        self._cache.dump('populations', new_populations)
        self._cache.dump('ref_populations', ref_populations)
        del self._local_cache
