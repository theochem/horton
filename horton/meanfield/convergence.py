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
'''Evaluation of convergence criteria

   These implementations are independent of the SCF algorithms and can be used
   to double check convergence.
'''


import numpy as np

from horton.meanfield.wfn import RestrictedWFN, UnrestrictedWFN


__all__ = [
    'convergence_error_eigen', 'compute_commutator',
    'convergence_error_commutator'
]


def convergence_error_eigen(ham, wfn, lf, overlap):
    '''Compute the self-consistency error

       **Arguments:**

       ham
            A Hamiltonian instance.

       wfn
            The wavefunction to be teste.

       lf
            The linalg factory to be used.

       overlap
            The overlap operator.

       **Returns:**

       error
            The SCF error. This measure (not this function) is also used
            in some SCF algorithms to check for convergence.
    '''
    if isinstance(wfn, RestrictedWFN):
        fock = lf.create_one_body()
        # Construct the Fock operator
        ham.compute_fock(fock, None)
        # Compute error
        return lf.error_eigen(fock, overlap, wfn.exp_alpha)
    elif isinstance(wfn, UnrestrictedWFN):
        fock_alpha = lf.create_one_body()
        fock_beta = lf.create_one_body()
        # Construct the Fock operators
        ham.compute_fock(fock_alpha, fock_beta)
        # Compute errors
        error_alpha = lf.error_eigen(fock_alpha, overlap, wfn.exp_alpha)
        error_beta = lf.error_eigen(fock_beta, overlap, wfn.exp_beta)
        return max(error_alpha, error_beta)
    else:
        raise NotImplementedError


def compute_commutator(dm, fock, overlap, work, output):
    '''Compute the dm-fock commutator, including an overlap matrix

       **Arguments:** (all OneBody objects)

       dm
            A density matrix

       fock
            A fock matrix

       overlap
            An overlap matrix

       work
            A temporary matrix

       output
            The output matrix in which the commutator, S.D.F-F.D.S, is stored.
    '''
    # construct sdf
    work.assign(overlap)
    work.idot(dm)
    work.idot(fock)
    output.assign(work)
    # construct fds and subtract
    work.assign(fock)
    work.idot(dm)
    work.idot(overlap)
    output.iadd(work, factor=-1)


def convergence_error_commutator(ham, wfn, lf, overlap):
    '''Compute the commutator error

       **Arguments:**

       ham
            A Hamiltonian instance.

       wfn
            The wavefunction to be teste.

       lf
            The linalg factory to be used.

       overlap
            The overlap operator.

       **Returns:**

       error
            The commutator error. This measure (not this function) is also used
            in some SCF algorithms to check for convergence.
    '''
    work = lf.create_one_body()
    if isinstance(wfn, RestrictedWFN):
        fock = lf.create_one_body()
        commutator = lf.create_one_body()
        # Construct the Fock operator
        ham.compute_fock(fock, None)
        # Compute commutator
        compute_commutator(wfn.dm_alpha, fock, overlap, work, commutator)
        # Compute norm
        normsq = commutator.expectation_value(commutator)
        return np.sqrt(normsq)
    elif isinstance(wfn, UnrestrictedWFN):
        fock_alpha = lf.create_one_body()
        fock_beta = lf.create_one_body()
        # Construct the Fock operators
        ham.compute_fock(fock_alpha, fock_beta)
        # Compute stuff for alpha
        compute_commutator(wfn.dm_alpha, fock_alpha, overlap, work, commutator)
        normsq_alpha = commutator.expectation_value(commutator)
        # Compute stuff for beta
        compute_commutator(wfn.dm_beta, fock_beta, overlap, work, commutator)
        normsq_beta = commutator.expectation_value(commutator)
        return np.sqrt(max(normsq_alpha, normsq_beta))
    else:
        raise NotImplementedError
