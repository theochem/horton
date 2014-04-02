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

from horton.log import log, timer
from horton.meanfield.scf_diis import DIISHistory, converge_scf_diis_cs
from horton.meanfield.wfn import RestrictedWFN
from horton.meanfield.hamiltonian import RestrictedEffectiveHamiltonian
from horton.quadprog import solve_safe


__all__ = ['converge_scf_cdiis']


@timer.with_section('SCF')
def converge_scf_cdiis(ham, wfn, lf, overlap, maxiter=128, threshold=1e-6, nvector=6, prune_old_states=False, skip_energy=False, scf_step='regular'):
    '''Minimize the energy of the wavefunction with the CDIIS algorithm

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       maxiter
            The maximum number of iterations. When set to None, the SCF loop
            will go one until convergence is reached.

       threshold
            The convergence threshold for the wavefunction

       prune_old_states
            When set to True, old states are pruned from the history when their
            coefficient is zero. Pruning starts at the oldest state and stops
            as soon as a state is encountered with a non-zero coefficient. Even
            if some newer states have a zero coefficient.

       skip_energy
            When set to True, the final energy is not computed.

       scf_step
            The type of SCF step to take after the interpolated states was
            create from the DIIS history. This can be 'regular', 'oda2' or
            'oda3'.

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.

       **Returns:** the number of iterations
    '''
    log.cite('pulay1980', 'using the commutator DIIS SCF algorithm')
    if isinstance(wfn, RestrictedWFN):
        assert isinstance(ham, RestrictedEffectiveHamiltonian)
        return converge_scf_cdiis_cs(ham, wfn, lf, overlap, maxiter, threshold, nvector, prune_old_states, skip_energy, scf_step)
    else:
        raise NotImplementedError


def converge_scf_cdiis_cs(ham, wfn, lf, overlap, maxiter=128, threshold=1e-6, nvector=6, prune_old_states=False, skip_energy=False, scf_step='regular'):
    '''Minimize the energy of the closed-shell wavefunction with CDIIS

       **Arguments:**

       ham
            A Hamiltonian instance.

       **Optional arguments:**

       maxiter
            The maximum number of iterations. When set to None, the SCF loop
            will go one until convergence is reached.

       threshold
            The convergence threshold for the wavefunction

       prune_old_states
            When set to True, old states are pruned from the history when their
            coefficient is zero. Pruning starts at the oldest state and stops
            as soon as a state is encountered with a non-zero coefficient. Even
            if some newer states have a zero coefficient.

       skip_energy
            When set to True, the final energy is not computed.

       scf_step
            The type of SCF step to take after the interpolated states was
            create from the DIIS history. This can be 'regular', 'oda2' or
            'oda3'.

       **Raises:**

       NoSCFConvergence
            if the convergence criteria are not met within the specified number
            of iterations.

       **Returns:** the number of iterations
    '''
    log.cite('pulay1980', 'the use of the commutator DIIS method')
    return converge_scf_diis_cs(ham, wfn, lf, overlap, PulayDIISHistory, maxiter, threshold, nvector, prune_old_states, skip_energy, scf_step)


class PulayDIISHistory(DIISHistory):
    '''A Pulay DIIS history object that keeps track of previous SCF solutions

       This type of DIIS is also called commutator DIIS, hence the name CDIIS.
    '''
    name = 'CDIIS'
    need_energy = False

    def __init__(self, lf, nvector, overlap):
        '''
           **Arguments:**

           lf
                The LinalgFactor used to create the one-body operators.

           nvector
                The maximum size of the history.

           overlap
                The overlap matrix.
        '''
        self.cdots = np.empty((nvector, nvector))
        self.cdots.fill(np.nan)
        DIISHistory.__init__(self, lf, nvector, overlap, [self.cdots])

    def _complete_cdots_matrix(self):
        '''Complete the matrix of dot products between commutators

           Even after multiple additions, this routine will fill up all the
           missing dot products in self.cdots.
        '''
        for i0 in xrange(self.nused-1, -1, -1):
            state0 = self.stack[i0]
            self.cdots[i0,i0] = state0.norm
            # Compute off-diagonal coefficients
            for i1 in xrange(i0):
                if np.isfinite(self.cdots[i0,i1]):
                    return
                state1 = self.stack[i1]
                cdot = state0.commutator.expectation_value(state1.commutator)
                self.cdots[i0,i1] = cdot
                self.cdots[i1,i0] = cdot

    def solve(self, dm_output, fock_output):
        '''Extrapolate a new density and/or fock matrix that should have the smallest commutator norm.

           **Arguments:**

           dm_output
                The output for the density matrix. If set to None, this is
                argument is ignored.

           fock_output
                The output for the Fock matrix. If set to None, this is
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
        self._build_combinations(coeffs, dm_output, fock_output)
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
