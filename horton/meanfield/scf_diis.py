# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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
"""Abstract DIIS code used by the different DIIS implementations"""

import numpy as np

from .exceptions import NoSCFConvergence
from .convergence import convergence_error_commutator
from .orbitals import Orbitals, check_dm
from .utils import compute_commutator

__all__ = []


class DIISSCFSolver(object):
    """Base class for all DIIS SCF solvers"""
    kind = 'dm'  # input/output variable is the density matrix

    def __init__(self, DIISHistoryClass, threshold=1e-6, maxiter=128, nvector=6, skip_energy=False,
                 prune_old_states=False):
        """
           **Arguments:**

           DIISHistoryClass
                A DIIS history class.

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
        """
        self.DIISHistoryClass = DIISHistoryClass
        self.threshold = threshold
        self.maxiter = maxiter
        self.nvector = nvector
        self.skip_energy = skip_energy
        self.prune_old_states = prune_old_states
        # Attributes where local variables of the __call__ method are stored:
        self._history = None
        self._focks = []
        self._orbs = []
        self.commutator = None
        self.biblio = []

    def __call__(self, ham, overlap, occ_model, *dms):
        """Find a self-consistent set of density matrices.

           **Arguments:**

           ham
                An effective Hamiltonian.

           overlap
                The overlap operator.

           occ_model
                Model for the orbital occupations.

           dm1, dm2, ...
                The initial density matrices. The number of dms must match
                ham.ndm.
        """
        # Some type checking
        if ham.ndm != len(dms):
            raise TypeError(
                'The number of initial density matrices does not match the Hamiltonian.')

        # Check input density matrices.
        for i in xrange(ham.ndm):
            check_dm(dms[i], overlap)
        occ_model.check_dms(overlap, *dms)

        # keep local variables as attributes for inspection/debugging by caller
        self._history = self.DIISHistoryClass(self.nvector, ham.ndm, ham.deriv_scale, overlap)
        self._focks = [np.zeros(overlap.shape) for i in xrange(ham.ndm)]
        self._orbs = [Orbitals(overlap.shape[0]) for i in xrange(ham.ndm)]

        print('5: Starting restricted closed-shell %s-SCF' % self._history.name)
        print("5: " + "-" * 70)
        print('5: Iter         Error        CN         Last nv Method          Energy       Change')
        print("5: " + "-" * 70)

        converged = False
        counter = 0
        while self.maxiter is None or counter < self.maxiter:
            # Construct the Fock operator from scratch if the history is empty:
            if self._history.nused == 0:
                # feed the latest density matrices in the hamiltonian
                ham.reset(*dms)
                # Construct the Fock operators
                ham.compute_fock(*self._focks)
                # Compute the energy if needed by the history
                energy = ham.compute_energy() if self._history.need_energy \
                    else None
                # Add the current fock+dm pair to the history
                error = self._history.add(energy, dms, self._focks)

                # Screen logging
                print('7:          DIIS add')
                if error < self.threshold:
                    converged = True
                    break
                print("7: ")
                energy_str = ' ' * 20 if energy is None else '% 20.13f' % energy
                print('5: %4i %12.5e                         %2i   %20s' % (
                    counter, error, self._history.nused, energy_str
                ))
                print("7: ")

            # Take a regular SCF step using the current fock matrix. Then
            # construct a new density matrix and fock matrix.
            for i in xrange(ham.ndm):
                self._orbs[i].from_fock(self._focks[i], overlap)
            occ_model.assign(*self._orbs)
            for i in xrange(ham.ndm):
                dms[i][:] = self._orbs[i].to_dm()
            ham.reset(*dms)
            energy = ham.compute_energy() if self._history.need_energy else None
            ham.compute_fock(*self._focks)

            # Add the current (dm, fock) pair to the history
            print('7:          DIIS add')
            error = self._history.add(energy, dms, self._focks)

            # break when converged
            if error < self.threshold:
                converged = True
                break

            # Screen logging
            print("7: ")
            energy_str = ' ' * 20 if energy is None else '% 20.13f' % energy
            print('5: %4i %12.5e                         %2i   %20s' % (
                counter, error, self._history.nused, energy_str
            ))
            print("7: ")

            # get extra/inter-polated Fock matrix
            while True:
                # The following method writes the interpolated dms and focks
                # in-place.
                energy_approx, coeffs, cn, method, error = self._history.solve(dms, self._focks)
                # if the error is small on the interpolated state, we have
                # converged to a solution that may have fractional occupation
                # numbers.
                if error < self.threshold:
                    converged = True
                    break
                # if coeffs[coeffs<0].sum() < -1:
                #    print('7:          DIIS (coeffs too negative) -> drop %i and retry' % self._history.stack[0].identity)
                #    self._history.shrink()
                if self._history.nused <= 2:
                    break
                if coeffs[-1] == 0.0:
                    print(
                        '7:          DIIS (last coeff zero) -> drop %i and retry' %
                        self._history.stack[
                            0].identity)
                    self._history.shrink()
                else:
                    break

            if False and len(coeffs) == 2:
                dms_tmp = [dm.copy() for dm in dms]
                xs = np.linspace(0.0, 1.0, 25)
                a, b = self._history._setup_equations()
                energies1 = []
                energies2 = []
                for x in xs:
                    x_coeffs = np.array([1 - x, x])
                    energies1.append(np.dot(x_coeffs, 0.5 * np.dot(a, x_coeffs) - b))
                    self._history._build_combinations(x_coeffs, dms_tmp, None)
                    ham.reset(*dms_tmp)
                    energies2.append(ham.compute_energy())
                    print(x, energies1[-1], energies2[-1])
                pt.clf()
                pt.plot(xs, energies1, label='est')
                pt.plot(xs, energies2, label='ref')
                pt.axvline(coeffs[1], color='k')
                pt.legend(loc=0)
                pt.savefig('diis_test_%05i.png' % counter)

            if energy_approx is not None:
                energy_change = energy_approx - min(state.energy for state in self._history.stack)
            else:
                energy_change = None

            # log
            self._history.log(coeffs)

            change_str = ' ' * 10 if energy_change is None else '% 12.7f' % energy_change
            print('5: %4i              %10.3e %12.7f %2i %s                      %12s' % (
                counter, cn, coeffs[-1], self._history.nused, method,
                change_str
            ))
            print("7: ")

            if self.prune_old_states:
                # get rid of old states with zero coeff
                for i in xrange(self._history.nused):
                    if coeffs[i] == 0.0:
                        print('7:          DIIS insignificant -> drop %i' % self._history.stack[
                            0].identity)
                        self._history.shrink()
                    else:
                        break

            # counter
            counter += 1

        if converged:
            print('5: %4i %12.5e (converged)' % (counter, error))
        print("7: ")

        if not self.skip_energy or self._history.need_energy:
            if not self._history.need_energy:
                ham.compute_energy()
            print(ham)

        if not converged:
            raise NoSCFConvergence

        return counter

    def error(self, ham, lf, overlap, *dms):
        return convergence_error_commutator(ham, lf, overlap, *dms)


class DIISState(object):
    """A single record (vector) in a DIIS history object."""

    def __init__(self, ndm, overlap):
        """
           **Arguments:**

           ndm
                The number of density matrices (and fock matrices) in one
                state.
           overlap
                The overlap matrix.
        """
        # Not all of these need to be used.
        self.ndm = ndm
        self.overlap = overlap
        self.energy = np.nan
        self.normsq = np.nan
        self.dms = [np.zeros(overlap.shape) for i in xrange(self.ndm)]
        self.focks = [np.zeros(overlap.shape) for i in xrange(self.ndm)]
        self.commutators = [np.zeros(overlap.shape) for i in xrange(self.ndm)]
        self.identity = None  # every state has a different id.

    def clear(self):
        """Reset this record."""
        self.energy = np.nan
        self.normsq = np.nan
        for i in xrange(self.ndm):
            self.dms[i][:] = 0.0
            self.focks[i][:] = 0.0
            self.commutators[i][:] = 0.0

    def assign(self, identity, energy, dms, focks):
        """Assign a new state.

           **Arguments:**

           identity
                A unique id for the new state.

           energy
                The energy of the new state.

           dm
                The density matrix of the new state.

           fock
                The Fock matrix of the new state.
        """
        self.identity = identity
        self.energy = energy
        self.normsq = 0.0
        for i in xrange(self.ndm):
            self.dms[i][:] = dms[i]
            self.focks[i][:] = focks[i]
            self.commutators[i][:] = compute_commutator(dms[i], focks[i], self.overlap)
            self.normsq += np.einsum('ab,ab', self.commutators[i], self.commutators[i])


class DIISHistory(object):
    """A base class of DIIS histories"""
    name = None
    need_energy = None

    def __init__(self, nvector, ndm, deriv_scale, overlap, dots_matrices):
        """
           **Arguments:**

           nvector
                The maximum size of the history.

           ndm
                The number of density matrices (and fock matrices) in one
                state.

           deriv_scale
                The deriv_scale attribute of the Effective Hamiltonian

           overlap
                The overlap matrix.

           dots_matrices
                Matrices in which dot products will be stored

           **Useful attributes:**

           used
                The actual number of vectors in the history.
        """
        self.stack = [DIISState(ndm, overlap) for i in xrange(nvector)]
        self.ndm = ndm
        self.deriv_scale = deriv_scale
        self.overlap = overlap
        self.dots_matrices = dots_matrices
        self.nused = 0
        self.idcounter = 0

    def _get_nvector(self):
        """The maximum size of the history"""
        return len(self.stack)

    nvector = property(_get_nvector)

    def log(self, coeffs):
        eref = min(state.energy for state in self.stack[:self.nused])
        if eref is None:
            print('5:          DIIS history          normsq       coeff         id')
            for i in xrange(self.nused):
                state = self.stack[i]
                print('5:          DIIS history  %12.5e  %12.7f   %8i' % (
                    state.normsq, coeffs[i], state.identity))
        else:
            print('5:          DIIS history          normsq      energy         coeff         id')
            for i in xrange(self.nused):
                state = self.stack[i]
                print('5:          DIIS history  %12.5e  %12.5e  %12.7f   %8i' % (
                    state.normsq, state.energy - eref, coeffs[i], state.identity))
        print("5: ")

    def solve(self, dms_output, focks_output):
        """Inter- or extrapolate new density and/or fock matrices.

           **Arguments:**

           dms_output
                The output for the density matrices. If set to None, this is
                argument is ignored.

           focks_output
                The output for the Fock matrices. If set to None, this is
                argument is ignored.
        """
        raise NotImplementedError

    def shrink(self):
        """Remove the oldest item from the history"""
        self.nused -= 1
        state = self.stack.pop(0)
        state.clear()
        self.stack.append(state)
        for dots in self.dots_matrices:
            dots[:-1] = dots[1:]
            dots[:, :-1] = dots[:, 1:]
            dots[-1] = np.nan
            dots[:, -1] = np.nan

    def add(self, energy, dms, focks):
        """Add new state to the history.

           **Arguments:**

           energy
                The energy of the new state.

           dms
                A list of density matrices of the new state.

           focks
                A list of Fock matrix of the new state.

           **Returns**: the square root of commutator error for the given pairs
           of density and Fock matrices.
        """
        if len(dms) != self.ndm or len(focks) != self.ndm:
            raise TypeError('The number of density and Fock matrices must match the ndm parameter.')
        # There must be a free spot. If needed, make one.
        if self.nused == self.nvector:
            self.shrink()

        # assign dm and fock
        state = self.stack[self.nused]
        state.assign(self.idcounter, energy, dms, focks)
        self.idcounter += 1

        # prepare for next iteration
        self.nused += 1
        return np.sqrt(state.normsq)

    def _build_combinations(self, coeffs, dms_output, focks_output):
        """Construct a linear combination of density/fock matrices

           **Arguments:**

           coeffs
                The linear mixing coefficients for the previous SCF states.

           dms_output
                A list of output density matrices. (Ignored if None)

           focks_output
                A list of output density matrices. (Ignored if None)

           **Returns:** the commutator error, only when both dms_output and
           focks_output are given.
        """
        if dms_output is not None:
            if len(dms_output) != self.ndm:
                raise TypeError('The number of density matrices must match the ndm parameter.')
            for i in xrange(self.ndm):
                dms_stack = [self.stack[j].dms[i] for j in xrange(self.nused)]
                self._linear_combination(coeffs, dms_stack, dms_output[i])
        if focks_output is not None:
            if len(focks_output) != self.ndm:
                raise TypeError('The number of Fock matrices must match the ndm parameter.')
            for i in xrange(self.ndm):
                focks_stack = [self.stack[j].focks[i] for j in xrange(self.nused)]
                self._linear_combination(coeffs, focks_stack, focks_output[i])
        if not (dms_output is None or focks_output is None):
            errorsq = 0.0
            for i in xrange(self.ndm):
                self.commutator = compute_commutator(dms_output[i], focks_output[i], self.overlap)
                errorsq += np.einsum('ab,ab', self.commutator, self.commutator)
            return errorsq ** 0.5

    def _linear_combination(self, coeffs, ops, output):
        """Make a linear combination of two-index objects

           **Arguments:**

           coeffs
                The linear mixing coefficients for the previous SCF states.

           ops
                A list of input operators.

           output
                The output operator.
        """
        output[:] = 0
        for i in xrange(self.nused):
            output += coeffs[i] * ops[i]
