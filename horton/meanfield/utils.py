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

__all__ = [
    'get_level_shift', 'get_spin', 'get_homo_lumo',
    'compute_commutator',
]

boltzmann = 3.1668154051341965e-06


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


def check_type(name, instance, *Classes):
    """Check type of argument with given name against list of types

       **Arguments:**

       name
            The name of the argument being checked.

       instance
            The object being checked.

       Classes
            A list of allowed types.
    """
    if len(Classes) == 0:
        raise TypeError('Type checking with an empty list of classes. This is a simple bug!')
    match = False
    for Class in Classes:
        if isinstance(instance, Class):
            match = True
            break
    if not match:
        classes_parts = ['\'', Classes[0].__name__, '\'']
        for Class in Classes[1:-1]:
            classes_parts.extend([', ``', Class.__name__, '\''])
        if len(Classes) > 1:
            classes_parts.extend(['or \'', Class.__name__, '\''])
        raise TypeError('The argument \'%s\' must be an instance of %s. Got a \'%s\' instance instead.' % (
            name, ''.join(classes_parts), instance.__class__.__name__
        ))


def fac2(n):
    result = 1
    while n > 1:
        result *= n
        n -= 2
    return result
