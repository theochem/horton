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
from horton.cache import just_once
from horton.log import log


__all__ = ['HirshfeldCPart']


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

    @just_once
    def _init_at_weights(self):
        # A) Compute all the pro-atomic densities.
        pro_atoms = []
        for i in xrange(self._system.natom):
            if log.do_medium:
                log('Computing pro-atom %i (%i)' % (i, self._system.natom))
            pro_atom = np.zeros(self._ui_grid.shape)
            log.mem.announce(pro_atom.nbytes)
            spline = self._proatomdb.get_hirshfeld_proatom_fn(self._system.numbers[i])
            center = self._system.coordinates[i]
            self._ui_grid.eval_spline(spline, center, pro_atom)
            pro_atom += 1e-100 # to avoid division by zero
            pro_atoms.append(pro_atom)
            del pro_atom

        # Compute the pro-molecule.
        pro_molecule = np.zeros(self._ui_grid.shape)
        log.mem.announce(pro_molecule.nbytes)
        for i in xrange(self._system.natom):
            pro_molecule += pro_atoms[i]

        # Construct the atomic weights
        for i in xrange(self._system.natom):
            pro_atoms[i] /= pro_molecule
            self._cache.dump('at_weights', i, pro_atoms[i], own=True)

        # Construct the relative molecular density
        pro_molecule *= -1
        pro_molecule += self._cache.load('mol_dens')
        self._cache.dump('rel_mol_dens', pro_molecule, own=True)

        # Store reference populations
        self._cache.dump('ref_populations', self._system.numbers.astype(float))
