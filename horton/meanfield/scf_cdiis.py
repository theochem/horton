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
'''Commutator DIIS Self-Consistent Field algorithm'''


import numpy as np

from horton.log import log
from horton.meanfield.scf_diis import DIISHistory, DIISSCFSolver
from horton.quadprog import solve_safe


__all__ = ['CDIISSCFSolver']


class CDIISSCFSolver(DIISSCFSolver):
    def __init__(self, threshold=1e-6, maxiter=128, nvector=6, skip_energy=False, prune_old_states=False):
        log.cite('pulay1980', 'the commutator DIIS SCF algorithm')
        DIISSCFSolver.__init__(self, CDIISHistory, threshold, maxiter, nvector, skip_energy, prune_old_states)


class CDIISHistory(DIISHistory):
    '''A commutator DIIS history object that keeps track of previous SCF solutions

       This type of DIIS is also called Pulay DIIS.
    '''
    name = 'CDIIS'
    need_energy = False

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
        self.cdots = np.empty((nvector, nvector))
        self.cdots.fill(np.nan)
        DIISHistory.__init__(self, lf, nvector, ndm, overlap, [self.cdots])

    def _complete_cdots_matrix(self):
        '''Complete the matrix of dot products between commutators

           Even after multiple additions, this routine will fill up all the
           missing dot products in self.cdots.
        '''
        for i0 in xrange(self.nused-1, -1, -1):
            state0 = self.stack[i0]
            self.cdots[i0,i0] = state0.normsq
            # Compute off-diagonal coefficients
            for i1 in xrange(i0):
                if np.isfinite(self.cdots[i0,i1]):
                    return
                state1 = self.stack[i1]
                cdot = 0.0
                for j in xrange(self.ndm):
                    cdot += state0.commutators[j].expectation_value(state1.commutators[j])
                self.cdots[i0,i1] = cdot
                self.cdots[i1,i0] = cdot

    def solve(self, dms_output, focks_output):
        '''Extrapolate a new density and/or fock matrices that should have the smallest commutator norm.

           **Arguments:**

           dms_output
                The output for the density matrices. If set to None, this is
                argument is ignored.

           focks_output
                The output for the Fock matrices. If set to None, this is
                argument is ignored.
        '''
        # extrapolation only makes sense if there are two points
        assert self.nused >= 2
        # Fill in the missing commutators
        self._complete_cdots_matrix()
        coeffs = solve_cdiis(self.cdots[:self.nused,:self.nused])
        # get a condition number
        absevals = abs(np.linalg.eigvalsh(self.cdots[:self.nused,:self.nused]))
        cn = absevals.max()/absevals.min()
        # assign extrapolated fock
        self._build_combinations(coeffs, dms_output, focks_output)
        return None, coeffs, cn, 'C'


def solve_cdiis(a):
    r'''Solve the linear equations found in the cdiis method

       The following is minimized:

       .. math:
            \frac{1}{2} x^T a x

        under the constraint :math:`\sum_i x_i = 1`.

       **Arguments:**

       a
            The matrix a, an array of size (N,N).
    '''
    n = len(a)
    assert a.shape == (n, n)
    assert (a == a.T).all()
    a2 = np.zeros((n+1, n+1))
    a2[:n,:n] = a
    a2[n,:n] = 1
    a2[:n,n] = 1
    b2 = np.zeros(n+1)
    b2[n] = 1
    x2 = solve_safe(a2, b2)
    return x2[:n]
