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
'''Energy DIIS Self-Consistent Field algorithm'''

import numpy as np

from horton.log import log, timer
from horton.exceptions import NoSCFConvergence
from horton.meanfield.scf_diis import DIISHistory, converge_scf_diis_cs
from horton.meanfield.wfn import RestrictedWFN
from horton.quadprog import QPSolver


__all__ = ['converge_scf_ediis']


@timer.with_section('SCF')
def converge_scf_ediis(ham, wfn, lf, overlap, maxiter=128, threshold=1e-6, nvector=6, prune_old_states=False, scf_step='regular'):
    '''Minimize the energy of the wavefunction with the EDIIS algorithm

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
    log.cite('kudin2002', 'using the energy DIIS SCF algorithm')
    if isinstance(wfn, RestrictedWFN):
        return converge_scf_ediis_cs(ham, wfn, lf, overlap, maxiter, threshold, nvector, prune_old_states, scf_step)
    else:
        raise NotImplementedError


def converge_scf_ediis_cs(ham, wfn, lf, overlap, maxiter=128, threshold=1e-6, nvector=6, prune_old_states=False, scf_step='regular'):
    '''Minimize the energy of the closed-shell wavefunction with EDIIS

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
    log.cite('kudin2002', 'the use of the EDIIS method.')
    return converge_scf_diis_cs(ham, wfn, lf, overlap, EnergyDIISHistory, maxiter, threshold, nvector, prune_old_states, scf_step)


class EnergyDIISHistory(DIISHistory):
    '''A Energy DIIS history object that keeps track of previous SCF solutions'''
    name = 'EDIIS'
    need_energy = True

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
        # A matrix with dot products of all density and fock matrices
        # Note that the dots matrix is not symmetric!
        self.edots = np.empty((nvector, nvector))
        self.edots.fill(np.nan)
        DIISHistory.__init__(self, lf, nvector, overlap, [self.edots])

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
                self.edots[i0,i1] = state0.fock.expectation_value(state1.dm)
                if i0 != i1:
                    # Note that this matrix is not symmetric!
                    self.edots[i1,i0] = state0.dm.expectation_value(state1.fock)

    def _setup_equations(self):
        b = np.zeros((self.nused, self.nused), float)
        e = np.zeros(self.nused, float)
        for i0 in xrange(self.nused):
            e[i0] = -self.stack[i0].energy
            for i1 in xrange(i0+1):
                b[i0, i1] = -(self.edots[i0,i0] + self.edots[i1,i1] - self.edots[i0,i1] - self.edots[i1,i0])
                if i0 != i1:
                    b[i1, i0] = b[i0, i1]
        return b, e

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
        # for debugging purposes
        try:
            qps.check_solution(coeffs)
        except:
            qps.log(guess)
            raise
        cn = qps.compute_cn(coeffs != 0.0)
        # assign extrapolated fock
        self._build_combinations(coeffs, dm_output, fock_output)
        return energy, coeffs, cn, 'E'
