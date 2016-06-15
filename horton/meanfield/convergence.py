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
'''Evaluation of convergence criteria

   These implementations are independent of the SCF algorithms and can be used
   to double check convergence.
'''


from horton.meanfield.utils import compute_commutator


__all__ = [
    'convergence_error_eigen', 'convergence_error_commutator',
]


def convergence_error_eigen(ham, lf, overlap, *exps):
    '''Compute the self-consistency error

       **Arguments:**

       ham
            A Hamiltonian instance.

       lf
            The linalg factory to be used.

       overlap
            The overlap operator.

       exp1, exp2, ...
            A list of wavefunction expansion objects. (The number must match
            ham.ndm.)

       **Returns:** The SCF error. This measure (not this function) is also used
       in some SCF algorithms to check for convergence.
    '''
    if len(exps) != ham.ndm:
        raise TypeError('Expecting %i expansions, got %i.' % (ham.ndm, len(exps)))
    dms = [exp.to_dm() for exp in exps]
    ham.reset(*dms)
    focks = [lf.create_two_index() for i in xrange(ham.ndm)]
    ham.compute_fock(*focks)
    error = 0.0
    for i in xrange(ham.ndm):
        error += exps[i].error_eigen(focks[i], overlap)
    return error


def convergence_error_commutator(ham, lf, overlap, *dms):
    '''Compute the commutator error

       **Arguments:**

       ham
            A Hamiltonian instance.

       lf
            The linalg factory to be used.

       overlap
            The overlap operator.

       dm1, dm2, ...
            A list of density matrices. The numbers of dms must match ham.ndm.

       **Returns:** The commutator error. This measure (not this function) is
       also used in some SCF algorithms to check for convergence.
    '''
    if len(dms) != ham.ndm:
        raise TypeError('Expecting %i density matrices, got %i.' % (ham.ndm, len(dms)))
    ham.reset(*dms)
    focks = [lf.create_two_index() for i in xrange(ham.ndm)]
    ham.compute_fock(*focks)
    error = 0.0
    work = lf.create_two_index()
    commutator = lf.create_two_index()
    errorsq = 0.0
    for i in xrange(ham.ndm):
        compute_commutator(dms[i], focks[i], overlap, work, commutator)
        errorsq += commutator.contract_two('ab,ab', commutator)
    return errorsq**0.5
