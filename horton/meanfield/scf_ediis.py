# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
'''Energy DIIS SCF algorithm'''

import numpy as np

from horton.log import log
from horton.exceptions import NoSCFConvergence
from horton.meanfield.scf_diis import DIISHistory, DIISSCFSolver
from horton.quadprog import QPSolver
from horton.utils import doc_inherit


__all__ = ['EDIISSCFSolver']


class EDIISSCFSolver(DIISSCFSolver):
    '''The Energy DIIS SCF solver [kudin2002]_'''

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
        log.cite('kudin2002', 'the EDIIS method.')
        DIISSCFSolver.__init__(self, EDIISHistory, threshold, maxiter, nvector, skip_energy, prune_old_states)


class EDIISHistory(DIISHistory):
    '''A Energy DIIS history object that keeps track of previous SCF solutions'''
    name = 'EDIIS'
    need_energy = True

    def __init__(self, lf, nvector, ndm, deriv_scale, overlap):
        '''
           **Arguments:**

           lf
                The LinalgFactor used to create the two-index operators.

           ndm
                The number of density matrices (and fock matrices) in one
                state.

           deriv_scale
                The deriv_scale attribute of the Effective Hamiltonian

           overlap
                The overlap matrix.
        '''
        # A matrix with dot products of all density and fock matrices
        # Note that the dots matrix is not symmetric!
        self.edots = np.empty((nvector, nvector))
        self.edots.fill(np.nan)
        DIISHistory.__init__(self, lf, nvector, ndm, deriv_scale, overlap, [self.edots])

    def _complete_edots_matrix(self):
        '''Complete the matrix of dot products between density and fock matrices

           Even after multiple additions, this routine will fill up all the
           missing dot products in self.edots.
        '''
        # This routine  even works after multiple additions.
        for i0 in xrange(self.nused-1, -1, -1):
            if np.isfinite(self.edots[i0,i0]):
                return
            # Compute off-diagonal coefficients
            state0 = self.stack[i0]
            for i1 in xrange(i0+1):
                state1 = self.stack[i1]
                self.edots[i0,i1] = 0.0
                for j in xrange(self.ndm):
                    self.edots[i0,i1] += state0.focks[j].contract_two('ab,ba', state1.dms[j])
                if i0 != i1:
                    # Note that this matrix is not symmetric!
                    self.edots[i1,i0] = 0.0
                    for j in xrange(self.ndm):
                        self.edots[i1,i0] += state1.focks[j].contract_two('ab,ba', state0.dms[j])

    def _setup_equations(self):
        '''Compute the equations for the quadratic programming problem.'''
        b = np.zeros((self.nused, self.nused), float)
        e = np.zeros(self.nused, float)
        for i0 in xrange(self.nused):
            e[i0] = -self.stack[i0].energy
            for i1 in xrange(i0+1):
                b[i0, i1] = -0.5*self.deriv_scale*(self.edots[i0,i0] + self.edots[i1,i1] - self.edots[i0,i1] - self.edots[i1,i0])
                if i0 != i1:
                    b[i1, i0] = b[i0, i1]
        return b, e

    @doc_inherit(DIISHistory)
    def solve(self, dms_output, focks_output):
        # interpolation only makes sense if there are two points
        assert self.nused >= 2
        # Fill in the missing commutators
        self._complete_edots_matrix()
        assert not np.isnan(self.edots[:self.nused,:self.nused]).any()
        # Setup the equations
        b, e = self._setup_equations()
        # Check if solving these equations makes sense.
        if b.max() - b.min() == 0 and e.max() - e.min() == 0:
            raise NoSCFConvergence('Convergence criteria too tight for EDIIS')
        # solve the quadratic programming problem
        qps = QPSolver(b, e, np.ones((1,self.nused)), np.array([1.0]), eps=1e-6)
        if self.nused < 10:
            energy, coeffs = qps.find_brute()
            guess = None
        else:
            guess = np.zeros(self.nused)
            guess[e.argmax()] = 1.0
            energy, coeffs = qps.find_local(guess, 1.0)
        # for debugging purposes (negligible computational overhead)
        try:
            qps.check_solution(coeffs)
        except:
            qps.log(guess)
            raise
        cn = qps.compute_cn(coeffs != 0.0)
        # assign extrapolated fock
        error = self._build_combinations(coeffs, dms_output, focks_output)
        return energy, coeffs, cn, 'E', error
