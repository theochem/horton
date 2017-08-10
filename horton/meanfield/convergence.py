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
"""Evaluation of convergence criteria

   These implementations are independent of the SCF algorithms and can be used
   to double check convergence.
"""

import numpy as np

from .utils import compute_commutator

__all__ = [
    'convergence_error_eigen', 'convergence_error_commutator',
]


def convergence_error_eigen(ham, overlap, *orbs):
    """Compute the self-consistency error

    Parameters
    ----------
    ham
        A Hamiltonian instance.
    overlap
        The overlap operator.
    orb1, orb2, ...
        A list of wavefunction expansion objects. (The number must match
        ham.ndm.)

    Returns
    -------
    error: float
        The SCF error. This measure (not this function) is also used in some SCF
        algorithms to check for convergence.
    """
    if len(orbs) != ham.ndm:
        raise TypeError('Expecting %i sets of orbitals, got %i.' % (ham.ndm, len(orbs)))
    dms = [orb.to_dm() for orb in orbs]
    ham.reset(*dms)
    focks = [np.zeros(dms[0].shape) for i in xrange(ham.ndm)]
    ham.compute_fock(*focks)
    error = 0.0
    for i in xrange(ham.ndm):
        error += orbs[i].error_eigen(focks[i], overlap)
    return error


def convergence_error_commutator(ham, overlap, *dms):
    """Compute the commutator error

    Parameters
    ----------
    ham
        A Hamiltonian instance.
    overlap
        The overlap operator.
    dm1, dm2, ...
        A list of density matrices. The numbers of dms must match ham.ndm.

    Returns
    -------
    error: float
        The commutator error. This measure (not this function) is also used in some SCF
        algorithms to check for convergence.
    """
    if len(dms) != ham.ndm:
        raise TypeError('Expecting %i density matrices, got %i.' % (ham.ndm, len(dms)))
    ham.reset(*dms)
    focks = [np.zeros(dms[0].shape) for i in xrange(ham.ndm)]
    ham.compute_fock(*focks)
    errorsq = 0.0
    for i in xrange(ham.ndm):
        commutator = compute_commutator(dms[i], focks[i], overlap)
        errorsq += np.einsum('ab,ab', commutator, commutator)
    return errorsq ** 0.5
