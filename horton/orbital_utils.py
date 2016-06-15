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
'''Utility functions for orbital modifications'''


from horton.log import timer
from horton.matrix.base import TwoIndex, Expansion
from horton.utils import check_type, check_options


__all__ = ['rotate_orbitals', 'compute_unitary_matrix', 'transform_integrals',
           'split_core_active']


def rotate_orbitals(*args):
    '''Rotate AO/MO (or MO/MO) coefficient matrix such that C` = C*U

       **Arguments:**

       args
            Rotation matrices (TwoIndex instance) and AO/MO coefficients
            (Expansion instance). rots and exps are ordered such that rot1
            corresponds to exp1, etc., i.e., rot1, rot2,..., exp1, exp2,...
    '''
    exps = []
    rotmat = []
    for arg in args:
        if isinstance(arg, TwoIndex):
            rotmat.append(arg)
        elif isinstance(arg, Expansion):
            exps.append(arg)
        else:
            raise TypeError('argument of unsupported type: %s' % arg)
    for i in xrange(len(exps)):
        exps[i].assign_dot(exps[i], rotmat[i])


def compute_unitary_matrix(kappa):
    '''Determine a unitary matrix from a skew-symmetric matrix K as
       U = exp(-K) by approximating U = 1 - K + 1/2 K^2 + O(3)

       **Arguments:**

       kappa
            A skew-symmetric matrix (TwoIndex instance)
    '''
    out = kappa.new()
    out.assign_diagonal(1.0)

    #
    # Approximate unitary matrix
    # U = exp(-K) by U = 1 - K + 1/2 K^2 + O(3)
    #
    out.iadd(kappa, -1.0)
    out.iadd_dot(kappa, kappa, 0.5)
    #
    # orthogonalization because approximate U matrix might not be unitary/orthogonal:
    #
    out.iortho()
    return out


@timer.with_section('Index Trans')
def transform_integrals(one, two, indextrans='tensordot', *exps):
    '''Update MO integrals. Returns list of transformed 1- and 2-electron
       integrals according to a list of expansion coefficients.

       **Arguments:**

       one
           One-electron integrals in the AO basis. A TwoIndex instance.

       two
           Two-electron integrals in the AO basis.

       **Optional arguments:**

       indextrans
           Choice of 4-index transformation. Default 'tensordot'.

       args
           The expansion coefficients.

    '''
    exp = []
    two_mo = []
    one_mo = []
    if not all([isinstance(i, Expansion) for i in exps]):
        raise TypeError('argument of unsupported type: %s' %i)
    for arg in exps:
        exp.append(arg)

    #
    # Loop over all possible AO/MO coefficients. We distinguish between
    #   * restricted orbitals: only one expansion instance (alpha=beta), returns
    #                          one set of 1- and 2-body integrals
    #   * unrestricted orbitals: two expansion instances (alpha, beta), returns
    #                            two sets of 1-body, and three sets of 2-body
    #                            integrals (one_alpha,alpha and one_beta,beta)
    #                            and (<alpha alpha|alpha alpha>, )
    #                            (<alpha beta|alpha beta>, <beta beta|beta beta>)
    #                            Note that order of alpha and beta is determined
    #                            by the order of the expansion instances
    #
    for i in xrange(len(exps)):
        for j in xrange(i, len(exps)):
            #
            # Transform 2-electron part
            # Note that integrals are stored using physics' convention <12|12>
            #
            out4ind = two.new()
            out4ind.assign_four_index_transform(two, exp[i], exp[j], exp[i], exp[j], indextrans)
            two_mo.append(out4ind)

        #
        # Transform 1-electron part
        #
        out2ind = one.new()
        out2ind.assign_two_index_transform(one, exps[i])
        one_mo.append(out2ind)
    return one_mo, two_mo


def split_core_active(one, two, ecore, orb, ncore, nactive, indextrans='tensordot'):
    '''Reduce a Hamiltonian to an active space

       Works only for restricted wavefunctions.

       **Arguments:**

       one/two
            One and two-electron integrals.

       ecore
            The core energy of the given Hamiltonian. In the case of a standard
            molecular system, this is the nuclear nuclear repulsion.

       orb
            The MO expansion coefficients. An Expansion instance. If None,
            integrals are assued to be already transformed into the mo basis
            and no transformation is carried out in this function.

       ncore
            The number of frozen core orbitals (int)

       nactive
            The number of active orbitals (int)

       **Optional arguments:**

       indextrans
            4-index transformation (str). One of ``tensordot``, ``einsum``

       **Returns** a tuple with three values:

       one_small
            The one-body operator in the small space

       two_small
            The two-body operator in the small space

       ecore
            The core energy, i.e. the sum of the given core energy and HF
            contributions from the core orbitals.
    '''
    #
    # Check type/option of arguments
    #
    check_type('ncore', ncore, int)
    check_type('nactive', nactive, int)
    check_options('indextrans', indextrans, 'tensordot', 'einsum')
    if ncore <= 0 or nactive <= 0:
        raise ValueError('ncore and nactive must be strictly positive.')
    if nactive+ncore > one.nbasis:
        raise ValueError('More active orbitals than basis functions.')

    #
    # Optional transformation to mo basis
    #
    if orb is None:
        one_mo = one
        two_mo = two
    else:
        # No need to check orb. This is done in transform_integrals function
        (one_mo,), (two_mo,) = transform_integrals(one, two, indextrans, orb)

    # Core energy
    norb = one.nbasis
    #   One body term
    ecore += 2*one_mo.trace(0, ncore, 0, ncore)
    #   Direct part
    ecore += two_mo.slice_to_two('abab->ab', None, 2.0, True, 0, ncore, 0, ncore, 0, ncore, 0, ncore).sum()
    #   Exchange part
    ecore += two_mo.slice_to_two('abba->ab', None,-1.0, True, 0, ncore, 0, ncore, 0, ncore, 0, ncore).sum()

    # Active space one-body integrals
    one_mo_corr = one_mo.new()
    #   Direct part
    two_mo.contract_to_two('abcb->ac', one_mo_corr, 2.0, True, 0, norb, 0, ncore, 0, norb, 0, ncore)
    #   Exchange part
    two_mo.contract_to_two('abbc->ac', one_mo_corr,-1.0, False, 0, norb, 0, ncore, 0, ncore, 0, norb)
    one_mo.iadd(one_mo_corr, 1.0)
    #one_mo.iadd_t(one_mo_corr, 1.0)

    #
    # Store in smaller n-index objects
    #
    one_mo_small = one_mo.copy(ncore, ncore+nactive, ncore, ncore+nactive)
    two_mo_small = two_mo.copy(ncore, ncore+nactive, ncore, ncore+nactive,
                               ncore, ncore+nactive, ncore, ncore+nactive)

    # Done
    return one_mo_small, two_mo_small, ecore
