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
'''EDIIS+DIIS Self-Consistent Field algorithm'''


import numpy as np

from horton.log import log
from horton.meanfield.scf_diis import DIISHistory, DIISSCFSolver
from horton.meanfield.scf_cdiis import CDIISHistory
from horton.meanfield.scf_ediis import EDIISHistory


__all__ = ['EDIIS2SCFSolver']


class EDIIS2SCFSolver(DIISSCFSolver):
    def __init__(self, threshold=1e-6, maxiter=128, nvector=6, skip_energy=False, prune_old_states=False):
        log.cite('kudin2002', 'the EDIIS method.')
        DIISSCFSolver.__init__(self, EDIIS2History, threshold, maxiter, nvector, skip_energy, prune_old_states)


class EDIIS2History(EDIISHistory, CDIISHistory):
    '''A EDIIS+DIIS history object that keeps track of previous SCF solutions

       This method uses EDIIS for the first iterations and switches to CDIIS
       to as soon as some initial degree of convergence is achieved.
    '''
    name = 'EDIIS+DIIS'
    need_energy = True

    def __init__(self, lf, nvector, ndm, overlap):
        '''
           **Arguments:**

           lf
                The LinalgFactor used to create the one-body operators.

           nvector
                The maximum size of the history.

           ndm
                The number of density matrices (and fock matrices) in one
                state.

           overlap
                The overlap matrix.
        '''
        # for the EDIIS part
        self.edots = np.empty((nvector, nvector))
        self.edots.fill(np.nan)
        # for the CDIIS part
        self.cdots = np.empty((nvector, nvector))
        self.cdots.fill(np.nan)
        DIISHistory.__init__(self, lf, nvector, ndm, overlap, [self.edots, self.cdots])

    def solve(self, dms_output, focks_output):
        '''Extrapolate a new density and/or fock matrix that should have the smallest commutator norm.

           **Arguments:**

           dms_output
                The output for the density matrices. If set to None, this is
                argument is ignored.

           focks_output
                The output for the Fock matrices. If set to None, this is
                argument is ignored.
        '''
        errmax = max(state.normsq for state in self.stack)
        if errmax > 1e-1:
            return EDIISHistory.solve(self, dms_output, focks_output)
        elif errmax < 1e-4:
            return CDIISHistory.solve(self, dms_output, focks_output)
        else:
            energy1, coeffs1, cn1, method1, error = CDIISHistory.solve(self, None, None)
            energy2, coeffs2, cn2, method2, error = EDIISHistory.solve(self, None, None)
            coeffs = 10*errmax*coeffs2 + (1-10*errmax)*coeffs1
            error = self._build_combinations(coeffs, dms_output, focks_output)
            return None, coeffs, max(cn1, cn2), 'M', error
