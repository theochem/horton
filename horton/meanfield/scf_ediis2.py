# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
'''EDIIS+DIIS SCF algorithm'''


import numpy as np

from horton.log import biblio
from horton.meanfield.scf_diis import DIISHistory, DIISSCFSolver
from horton.meanfield.scf_cdiis import CDIISHistory
from horton.meanfield.scf_ediis import EDIISHistory
from horton.utils import doc_inherit


__all__ = ['EDIIS2SCFSolver']


class EDIIS2SCFSolver(DIISSCFSolver):
    '''The EDIIS+DIIS SCF solver [kudin2002]_'''

    def __init__(self, threshold=1e-6, maxiter=128, nvector=6, skip_energy=False, prune_old_states=False):
        '''
           **Optional arguments:**

           maxiter
                The maximum number of iterations. When set to None, the SCF loop
                will go one until convergence is reached.

           threshold
                The convergence threshold for the wavefunction

           skip_energy
                When set to True, the final energy is not computed. Note that some
                DIIS variants need to compute the energy anyway. for these methods
                this option is irrelevant.

           prune_old_states
                When set to True, old states are pruned from the history when their
                coefficient is zero. Pruning starts at the oldest state and stops
                as soon as a state is encountered with a non-zero coefficient. Even
                if some newer states have a zero coefficient.
        '''
        biblio.cite('kudin2002', 'the EDIIS method.')
        DIISSCFSolver.__init__(self, EDIIS2History, threshold, maxiter, nvector, skip_energy, prune_old_states)


class EDIIS2History(EDIISHistory, CDIISHistory):
    '''A EDIIS+DIIS history object that keeps track of previous SCF solutions

       This method uses EDIIS for the first iterations and switches to CDIIS
       to as soon as some initial degree of convergence is achieved.
    '''
    name = 'EDIIS+DIIS'
    need_energy = True

    def __init__(self, lf, nvector, ndm, deriv_scale, overlap):
        '''
           **Arguments:**

           lf
                The LinalgFactor used to create the two-index operators.

           nvector
                The maximum size of the history.

           ndm
                The number of density matrices (and fock matrices) in one
                state.

           deriv_scale
                The deriv_scale attribute of the Effective Hamiltonian

           overlap
                The overlap matrix.
        '''
        # for the EDIIS part
        self.edots = np.empty((nvector, nvector))
        self.edots.fill(np.nan)
        # for the CDIIS part
        self.cdots = np.empty((nvector, nvector))
        self.cdots.fill(np.nan)
        DIISHistory.__init__(self, lf, nvector, ndm, deriv_scale, overlap, [self.edots, self.cdots])

    @doc_inherit(DIISHistory)
    def solve(self, dms_output, focks_output):
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
