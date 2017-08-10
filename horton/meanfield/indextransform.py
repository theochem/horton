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
"""Utility functions for orbital modifications."""


import numpy as np

__all__ = ['four_index_transform', 'transform_integrals', 'split_core_active',
           'four_index_transform_cholesky', 'transform_integrals_cholesky', 'split_core_active_cholesky']


def _parse_four_index_transform_orbs(orb0, orb1, orb2, orb3):
    """Parse the optional arguments exp1, exp2 and exp3.

    Parameters
    ----------
    orb0, orb1, orb2, orb3
        Four sets of orbitals for the mo transformation. Some may be None but only the
        following not None combinations are allowed:

        * ``(orb0,)``: maintain eight-fold symmetry (if any)
        * ``(orb0, orb1)``: maintain four-fold symmetry (if any)
        * ``(orb0, orb2)``: maintain two-fold symmetry (if any)
        * ``(orb0, orb1, orb2, orb3)``: break all symmetry

    Returns
    -------
    orb0, orb1, orb2, orb3. (All not None)
    """
    # Four supported situations
    if orb1 is None and orb2 is None and orb3 is None:
        # maintains eight-fold symmetry
        orb1 = orb0
        orb2 = orb0
        orb3 = orb0
    elif orb2 is None and orb3 is None:
        # maintains four-fold symmetry
        orb2 = orb0
        orb3 = orb1
    elif orb1 is None and orb3 is None:
        # maintains two-fold symmetry
        orb1 = orb0
        orb3 = orb2
    elif orb1 is None or orb2 is None or orb3 is None:
        # the only other allowed case is no symmetry.
        raise TypeError('It is not clear how to interpret the optional arguments orb1, orb2 and orb3.')
    return orb0, orb1, orb2, orb3


def four_index_transform(ao_integrals, orb0, orb1=None, orb2=None, orb3=None, method='tensordot'):
    """Perform four index transformation.

    Parameters
    ----------
    oa_integrals
        A four-index array with integrals over atomic orbitals.
    orb0
        A Orbitalas object with molecular orbitals
    orb1, orb2, orb3
        Can be provided to transform each index differently.
    method
        Either ``einsum`` or ``tensordot`` (default).

    Returns
    -------
    mo_integrals
        A four-index array with the integrals in the MO basis.
    """
    # parse arguments
    orb0, orb1, orb2, orb3 = _parse_four_index_transform_orbs(orb0, orb1, orb2, orb3)
    # actual transform
    result = np.zeros(ao_integrals.shape)
    if method == 'einsum':
        # The order of the dot products is according to literature
        # conventions.
        result[:] = np.einsum('sd,pqrs->pqrd', orb3.coeffs, ao_integrals, casting='no', order='C')
        result[:] = np.einsum('rc,pqrd->pqcd', orb2.coeffs, result, casting='no', order='C')
        result[:] = np.einsum('qb,pqcd->pbcd', orb1.coeffs, result, casting='no', order='C')
        result[:] = np.einsum('pa,pbcd->abcd', orb0.coeffs, result, casting='no', order='C')
    elif method == 'tensordot':
        # because the way tensordot works, the order of the dot products is
        # not according to literature conventions.
        result[:] = np.tensordot(ao_integrals, orb0.coeffs, axes=([0],[0]))
        result[:] = np.tensordot(result, orb1.coeffs, axes=([0],[0]))
        result[:] = np.tensordot(result, orb2.coeffs, axes=([0],[0]))
        result[:] = np.tensordot(result, orb3.coeffs, axes=([0],[0]))
    else:
        raise ValueError('The method must either be \'einsum\' or \'tensordot\'.')
    # Symmetrize the result
    result[:] = result + result.transpose(1,0,3,2)
    result[:] = result + result.transpose(2,3,0,1)
    result[:] = result + result.transpose(0,3,2,1)
    result /= 8
    return result


def transform_integrals(one, two, method='tensordot', *orbs):
    """Transform integrals to MO basis.

    Parameters
    ----------
    one
       One-electron integrals in the AO basis.
    two
       Two-electron integrals in the AO basis.
    method
       Choice of 4-index transformation: 'tensordot' or 'einsum'.
    orbs
       A list of Orbitals objects.

    Returns
    -------
    list of transformed 1- and 2-electron integrals according to a list of Orbitals
    instances. We distinguish between:

    * restricted orbitals: only one Ortbials instance (alpha=beta), returns
      one set of 1- and 2-body integrals
    * unrestricted orbitals: two Ortbials instances (alpha, beta), returns
      two sets of 1-body, and three sets of 2-body integrals:

      - one_{alpha,alpha} and
      - one_{beta,beta}
      - two_{alpha alpha|alpha alpha}
      - two_{<alpha beta|alpha beta}
      - two_{beta beta|beta beta}

      Note that order of alpha and beta is determined by the order of the Orbitals
      instances.
    """
    two_mo = []
    one_mo = []
    for i0, orb0 in enumerate(orbs):
        for orb1 in orbs[i0:]:
            # Note that integrals are stored using physics' convention <12|12>
            two_mo.append(four_index_transform(two, orb0, orb1, orb0, orb1, method))
        # Transform 1-electron part
        one_mo.append(reduce(np.dot, [orb0.coeffs.T, one, orb0.coeffs]))
    return one_mo, two_mo


def split_core_active(one, two, ecore, orb, ncore, nactive, indextrans='tensordot'):
    """Reduce a Hamiltonian to an active space.

    Works only for restricted wavefunctions.

    Parameters
    ----------
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
    indextrans
        4-index transformation (str). One of ``tensordot``, ``einsum``

    Returns
    -------
    one_small
        The one-body operator in the small space.
    two_small
        The two-body operator in the small space.
    ecore
        The core energy, i.e. the sum of the given core energy and HF
        contributions from the core orbitals.
    """
    # Check type/option of arguments
    if ncore < 0:
        raise ValueError('ncore must be positive.')
    if nactive <= 0:
        raise ValueError('ncore must be strictly positive.')
    if nactive + ncore > one.shape[0]:
        raise ValueError('More active orbitals than basis functions.')

    # Optional transformation to mo basis
    if orb is None:
        one_mo = one
        two_mo = two
    else:
        (one_mo,), (two_mo,) = transform_integrals(one, two, indextrans, orb)

    # Core energy
    #   One body term
    ecore += 2*np.trace(one_mo[:ncore, :ncore])
    #   Direct part
    ecore += 2*np.einsum('abab', two_mo[:ncore, :ncore, :ncore, :ncore])
    #   Exchange part
    ecore -= np.einsum('abba', two_mo[:ncore, :ncore, :ncore, :ncore])

    # Active space one-body integrals
    norb = ncore + nactive
    one_mo_small = one_mo[ncore:norb, ncore:norb].copy()
    #   Direct part
    one_mo_small += 2*np.einsum('abcb->ac', two_mo[ncore:norb, :ncore, ncore:norb, :ncore])
    #   Exchange part
    one_mo_small -= np.einsum('abbc->ac', two_mo[ncore:norb, :ncore, :ncore, ncore:norb])

    # Active space two-body integrals
    two_mo_small = two_mo[ncore:norb, ncore:norb, ncore:norb, ncore:norb]

    return one_mo_small, two_mo_small, ecore


def four_index_transform_cholesky(ao_integrals, orb0, orb1=None, method='tensordot'):
    """Perform four index transformation on a Cholesky-decomposed four-index object.

    Parameters
    ----------
    oa_integrals : np.ndarray, shape=(nvec, nbasis, nbasis)
        Cholesky decomposition of four-index object in the AO basis.
    orb0
        A Orbitals object with molecular orbitals.
    orb1
        Can be provided to transform the second index differently.
    method
        Either ``einsum`` or ``tensordot`` (default).
    """
    if orb1 is None:
        orb1 = orb0
    result = np.zeros(ao_integrals.shape)
    if method == 'einsum':
        result = np.einsum('ai,kac->kic', orb0.coeffs, ao_integrals)
        result = np.einsum('cj,kic->kij', orb1.coeffs, result)
    elif method == 'tensordot':
        result = np.tensordot(ao_integrals, orb0.coeffs, axes=([1],[0]))
        result = np.tensordot(result, orb1.coeffs, axes=([1],[0]))
    else:
        raise ValueError('The method must either be \'einsum\' or \'tensordot\'.')
    return result


def transform_integrals_cholesky(one, two, method='tensordot', *orbs):
    """Transform integrals to MO basis.

    Parameters
    ----------
    one
       One-electron integrals in the AO basis.
    two
       Cholesky decomposition of two-electron integrals in the AO basis.
    method
       Choice of 4-index transformation: 'tensordot' or 'einsum'.
    orbs
       A list of Orbitals objects.

    Returns
    -------
    list of transformed 1- and 2-electron integrals according to a list of Orbitals
    instances.  We distinguish between:

    * restricted orbitals: only one Ortbials instance (alpha=beta), returns
      one set of 1- and 2-body integrals
    * unrestricted orbitals: two expansion instances (alpha, beta), returns
      two sets of 1-body, and two sets of 2-body integrals.
    """
    two_mo = []
    one_mo = []
    for orb0 in orbs:
        two_mo.append(four_index_transform_cholesky(two, orb0, orb0, method))
        one_mo.append(reduce(np.dot, [orb0.coeffs.T, one, orb0.coeffs]))
    return one_mo, two_mo


def split_core_active_cholesky(one, two, ecore, orb, ncore, nactive, indextrans='tensordot'):
    """Reduce a Hamiltonian to an active space.

    Works only for restricted wavefunctions.

    Parameters
    ----------
    one/two
        One and two-electron integrals. A Cholesky decomposition of the two-electron
        integrals must be provided.
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
    indextrans
        4-index transformation (str). One of ``tensordot``, ``einsum``

    Returns
    -------
    one_small
        The one-body operator in the small space.
    two_small
        A Cholesky decomposition of the two-body operator in the small space.
    ecore
        The core energy, i.e. the sum of the given core energy and HF
        contributions from the core orbitals.
    """
    # Check type/option of arguments
    if ncore < 0:
        raise ValueError('ncore must be positive.')
    if nactive <= 0:
        raise ValueError('ncore must be strictly positive.')
    if nactive + ncore > one.shape[0]:
        raise ValueError('More active orbitals than basis functions.')

    # Optional transformation to mo basis
    if orb is None:
        one_mo = one
        two_mo = two
    else:
        (one_mo,), (two_mo,) = transform_integrals_cholesky(one, two, indextrans, orb)

    # Core energy
    #   One body term
    ecore += 2*np.trace(one_mo[:ncore, :ncore])
    #   Direct part
    ecore += 2*np.einsum('xaa,xbb', two_mo[:, :ncore, :ncore], two_mo[:, :ncore, :ncore])
    #   Exchange part
    ecore -= np.einsum('xab,xba', two_mo[:, :ncore, :ncore], two_mo[:, :ncore, :ncore])

    # Active space one-body integrals
    norb = ncore + nactive
    one_mo_small = one_mo[ncore:norb, ncore:norb].copy()
    #   Direct part
    one_mo_small += 2*np.einsum('xac,xbb->ac', two_mo[:, ncore:norb, ncore:norb], two_mo[:, :ncore, :ncore])
    #   Exchange part
    one_mo_small -= np.einsum('xab,xbc->ac', two_mo[:, ncore:norb, :ncore], two_mo[:, :ncore, ncore:norb])

    # Active space two-body integrals
    two_mo_small = two_mo[:, ncore:norb, ncore:norb]

    return one_mo_small, two_mo_small, ecore
