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
"""Basic Self-Consistent Field (SCF) algorithm."""

import numpy as np

from horton.exceptions import NoSCFConvergence
from .convergence import convergence_error_eigen
from .utils import get_level_shift

__all__ = ['PlainSCFSolver']


class PlainSCFSolver(object):
    """A bare-bones SCF solver without mixing."""

    kind = 'orb'  # input/output variable is the Orbitals object

    def __init__(self, threshold=1e-8, maxiter=128, skip_energy=False, level_shift=0.0):
        """Initialize a PlainSCFSolver.

        Parameters
        ----------

        threshold : float
                    The convergence threshold for the wavefunction.
        maxiter : int
                  The maximum number of iterations. When set to None, the SCF loop will go
                  one until convergence is reached.
        skip_energy : bool
                      When set to True, the final energy is not computed.
        level_shift : float
                      When set to non-zero, level-shifting is applied and the value of the
                      argument controls the magnitude of the level shift. This argument
                      cannot be negative.
        """
        self.maxiter = maxiter
        self.threshold = threshold
        self.skip_energy = skip_energy
        if level_shift < 0:
            raise ValueError('The level_shift argument cannot be negative.')
        self.level_shift = level_shift

    def __call__(self, ham, overlap, occ_model, *orbs):
        """Find a self-consistent set of orbitals.

        Parameters
        ----------

        ham : EffHam
            An effective Hamiltonian.
        overlap : np.ndarray, shape=(nbasis, nbasis)
            The overlap operator.
        occ_model : OccModel
            Model for the orbital occupations.
        orb1, orb2, ... : Orbitals
            The initial orbitals. The number of dms must match ham.ndm.
        """
        # Some type checking
        if ham.ndm != len(orbs):
            raise TypeError(
                'The number of initial orbital expansions does not match the Hamiltonian.')
        # Impose the requested occupation numbers
        occ_model.assign(*orbs)
        # Check the orthogonality of the orbitals
        for orb in orbs:
            orb.check_normalization(overlap)

        print('Starting plain SCF solver. ndm=%i' % ham.ndm)
        print("5: " + "-" * 70)
        print('Iter         Error')
        print("5: " + "-" * 70)

        focks = [np.zeros(overlap.shape) for i in xrange(ham.ndm)]
        dms = [None] * ham.ndm
        converged = False
        counter = 0
        while self.maxiter is None or counter < self.maxiter:
            # convert the orbital expansions to density matrices
            for i in xrange(ham.ndm):
                dms[i] = orbs[i].to_dm()
            # feed the latest density matrices in the hamiltonian
            ham.reset(*dms)
            # Construct the Fock operator
            ham.compute_fock(*focks)
            # Check for convergence
            error = 0.0
            for i in xrange(ham.ndm):
                error += orbs[i].error_eigen(focks[i], overlap)
            print('5: %4i  %12.5e' % (counter, error))
            if error < self.threshold:
                converged = True
                break
            # If requested, add the level shift to the Fock operator
            if self.level_shift > 0:
                for i in xrange(ham.ndm):
                    # The normal behavior is to shift down the occupied levels.
                    focks[i] += -self.level_shift * get_level_shift(dms[i], overlap)
            # Diagonalize the fock operators to obtain new orbitals and
            for i in xrange(ham.ndm):
                orbs[i].from_fock(focks[i], overlap)
                # If requested, compensate for level-shift. This compensation
                # is only correct when the SCF has converged.
                if self.level_shift > 0:
                    orbs[i].energies[:] += self.level_shift * orbs[i].occupations
            # Assign new occupation numbers.
            occ_model.assign(*orbs)
            # counter
            counter += 1

        print("5: ")

        if not self.skip_energy:
            ham.compute_energy()
            print(ham)

        if not converged:
            raise NoSCFConvergence

        return counter

    def error(self, ham, overlap, *orbs):
        """See :py:func:`horton.meanfield.convergence.convergence_error_eigen`."""
        return convergence_error_eigen(ham, overlap, *orbs)
