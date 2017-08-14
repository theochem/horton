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
"""Utility functions"""

import numpy as np

from .orbitals import Orbitals

__all__ = [
    'check_dm', 'get_level_shift', 'get_spin', 'get_homo_lumo',
    'compute_commutator',
]

boltzmann = 3.1668154051341965e-06


def check_dm(dm, overlap, eps=1e-4, occ_max=1.0):
    """Check if the density matrix has eigenvalues in the proper range.

    Parameters
    ----------
    dm : np.ndarray, shape=(nbasis, nbasis), dtype=float
        The density matrix
    overlap : np.ndarray, shape=(nbasis, nbasis), dtype=float
        The overlap matrix
    eps : float
        The threshold on the eigenvalue inequalities.
    occ_max : float
        The maximum occupation.

    Raises
    ------
    ValueError
        When the density matrix has wrong eigenvalues.
    """
    # construct natural orbitals
    orb = Orbitals(dm.shape[0])
    orb.derive_naturals(dm, overlap)
    if orb.occupations.min() < -eps:
        raise ValueError('The density matrix has eigenvalues considerably smaller than '
                         'zero. error=%e' % (orb.occupations.min()))
    if orb.occupations.max() > occ_max + eps:
        raise ValueError('The density matrix has eigenvalues considerably larger than '
                         'max. error=%e' % (orb.occupations.max() - 1))


def get_level_shift(dm, overlap):
    """Construct a level shift operator.

       **Arguments:**

       dm
            A density matrix.

       overlap
            The overlap matrix

       **Returns:** The level-shift operator.
    """
    return np.dot(overlap.T, np.dot(dm, overlap))


def get_spin(orb_alpha, orb_beta, overlap):
    """Returns the expectation values of the projected and squared spin

       **Arguments:**

       orb_alpha, orb_beta
            The alpha and beta orbitals.

       overlap
            The overlap matrix

       **Returns:** sz, ssq
    """
    nalpha = orb_alpha.occupations.sum()
    nbeta = orb_beta.occupations.sum()
    sz = (nalpha - nbeta) / 2
    correction = 0.0
    for ialpha in xrange(orb_alpha.nfn):
        if orb_alpha.occupations[ialpha] == 0.0:
            continue
        for ibeta in xrange(orb_beta.nfn):
            if orb_beta.occupations[ibeta] == 0.0:
                continue
            correction += np.dot(
                orb_alpha.coeffs[:, ialpha],
                np.dot(overlap, orb_beta.coeffs[:, ibeta])) ** 2

    ssq = sz * (sz + 1) + nbeta - correction
    print sz, ssq
    return sz, ssq


def get_homo_lumo(*orbs):
    """Return the HOMO and LUMO energy for the given expansion

       **Arguments:**

       orb1, orb2, ...
            Orbitals objects

       **Returns:** homo_energy, lumo_energy. (The second is None when all
       orbitals are occupied.)
    """
    homo_energy = max(orb.homo_energy for orb in orbs)
    lumo_energies = [orb.lumo_energy for orb in orbs]
    lumo_energies = [lumo_energy for lumo_energy in lumo_energies if lumo_energy is not None]
    if len(lumo_energies) == 0:
        lumo_energy = None
    else:
        lumo_energy = min(lumo_energies)
    return homo_energy, lumo_energy


def compute_commutator(dm, fock, overlap):
    """Compute the dm-fock commutator, including an overlap matrix

    Parameters
    ----------
    dm
        A density matrix
    fock
        A fock matrix
    overlap
        An overlap matrix

    Return
    ------
    commutator
    """
    return np.dot(overlap, np.dot(dm, fock)) - np.dot(fock, np.dot(dm, overlap))


def doc_inherit(base_class):
    """Docstring inheriting method descriptor

       doc_inherit decorator

       Usage:

       .. code-block:: python

            class Foo(object):
                def foo(self):
                    "Frobber"
                    pass

            class Bar(Foo):
                @doc_inherit(Foo)
                def foo(self):
                    pass

       Now, ``Bar.foo.__doc__ == Bar().foo.__doc__ == Foo.foo.__doc__ ==
       "Frobber"``
    """

    def decorator(method):
        overridden = getattr(base_class, method.__name__, None)
        if overridden is None:
            raise NameError('Can\'t find method \'%s\' in base class.')
        method.__doc__ = overridden.__doc__
        return method

    return decorator


# from horton.moments
def get_ncart_cumul(lmax):
    """The number of cartesian powers up to a given angular momentum, lmax."""
    return ((lmax+1)*(lmax+2)*(lmax+3))/6


# from horton.moments
def get_cartesian_powers(lmax):
    """Return an ordered list of power for x, y and z up to angular moment lmax

       **Arguments:**

       lmax
            The maximum angular momentum (0=s, 1=p, 2=d, ...)

       **Returns:** an array where each row corresponds to a multipole moment
       and each column corresponds to a power of x, y and z respectively. The
       rows are grouped per angular momentum, first s, them p, then d, and so
       on. Within one angular momentum the rows are sorted 'alphabetically',
       e.g. for l=2: xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz.
    """
    cartesian_powers = np.zeros((get_ncart_cumul(lmax), 3), dtype=int)
    counter = 0
    for l in xrange(0, lmax + 1):
        for nx in xrange(l + 1, -1, -1):
            for ny in xrange(l - nx, -1, -1):
                nz = l - ny - nx
                cartesian_powers[counter] = [nx, ny, nz]
                counter += 1
    return cartesian_powers


# from horton.moments
def rotate_cartesian_multipole(rmat, moments, mode):
    """Compute rotated Cartesian multipole moment/expansion.

       **Arguments:**

       rmat
            A (3,3) rotation matrix.

       moments
            A multipole moment/coeffs. The angular momentum is derived from the
            length of this vector.

       mode
            A string containing either 'moments' or 'coeffs'. In case if
            'moments', a Cartesian multipole moment rotation is carried out. In
            case of 'coeffs', the coefficients of a Cartesian multipole basis
            are rotated.

       **Returns:** rotated multipole.
    """
    l = ((9 + 8 * (len(moments) - 1)) ** 0.5 - 3) / 2
    if l - np.round(l) > 1e-10:
        raise ValueError('Could not determine l from number of moments.')
    l = int(np.round(l))

    if mode == 'coeffs':
        rcoeffs = rmat.T.ravel()
    elif mode == 'moments':
        rcoeffs = rmat.ravel()
    else:
        raise NotImplementedError
    result = np.zeros(len(moments))
    for i0 in xrange(len(moments)):
        rules = cartesian_transforms[l][i0]
        for rule in rules:
            i1 = rule[0]
            factor = rule[1]
            for j in rule[2:]:
                factor *= rcoeffs[j]
            if mode == 'coeffs':
                result[i1] += moments[i0] * factor
            elif mode == 'moments':
                result[i0] += moments[i1] * factor
            else:
                raise NotImplementedError
    return result
