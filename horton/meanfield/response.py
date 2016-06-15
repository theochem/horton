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
'''Evaluation of response functions'''

import numpy as np

from horton.log import timer


__all__ = ['compute_noninteracting_response']


@timer.with_section('KS Response')
def compute_noninteracting_response(exp, operators, work=None):
    '''Compute the non-interacting response matrix for a given orbital expansion

       **Arguments:**

       exp
            An instance of DenseExpansion.

       operators
            A list of one-body operators.

       **Optional arguments:**

       work
            A work array with shape (len(operators), nfn, nfn), where nfn is
            the number of occupied and virtual orbitals.

       **Returns:** a symmetric matrix where each element corresponds to a pair
       of operators. Note that this function is only reliable when no degenerate
       orbitals are present at the fermi level. For example, in case of
       fractional occupations in DFT, this method does not give the correct
       non-interacting response matrix.
    '''
    # Convert the operators to the orbital basis
    coeffs = exp.coeffs
    norb = exp.nfn
    nop = len(operators)

    if work is None:
        work = np.zeros((nop, norb, norb))
    for iop in xrange(nop):
        work[iop] = np.dot(coeffs.T, np.dot(operators[iop]._array, coeffs))

    # Put the orbital energies and the occupations in a convenient array
    energies = exp.energies
    occupations = exp.occupations
    with np.errstate(invalid='ignore'):
        prefacs = np.subtract.outer(occupations, occupations)/np.subtract.outer(energies, energies)
    # Purge divisions by zero. If degeneracies occur at the fermi level, this
    # way of computing the noninteracting response matrix is not correct anyway.
    #for iorb in xrange(norb):
    #    prefacs[iorb,iorb] = 0.0
    mask = occupations == occupations.reshape(-1, 1)
    mask |= energies == energies.reshape(-1, 1)
    prefacs[mask] = 0.0

    # Double loop over all pairs for operators. The diagonal element is computed
    # as a double check. Note that by construction the first term and its
    # complex conjugate of the SOS expressions,
    #
    #     X_s = \sum_\substrack{i \in \text{occ} \\ j \in \text{virt}}
    #               (< phi_i | A | phi_j > < phi_j | B | phi_i >)/(\epsilon_i - \epsilon_j)
    #           + c.c.
    #
    # are included. (The first term corresponds to the lower diagonal of prefacs
    # while the complex conjugate corresponds to the upper diagonal of prefacs.)
    result = np.zeros((nop, nop), float)
    for iop0 in xrange(nop):
        for iop1 in xrange(iop0+1):
            # evaluate the sum over states expression
            state_sum = (work[iop0]*work[iop1]*prefacs).sum()

            # store the state sum
            result[iop0, iop1] = state_sum
            result[iop1, iop0] = state_sum

    return result
