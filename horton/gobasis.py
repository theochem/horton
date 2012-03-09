# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>, ...
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
"""Molecular gaussian orbital basis set module.

   This module focuses on basis sets applied to a certain molecular structure,
   which is different from the typical database of atomic basis sets that are
   used to compose molecular basis sets.
"""


import numpy as np


__all__ = ['get_con_nbasis', 'GaussianOrbBasis']


def get_con_nbasis(con_type):
    """Return the number basis function for a given contraction type."""
    if con_type > 0:
        return (con_type+1)*(con_type+2)/2
    elif con_type == -1:
        raise ValueError
    else:
        return -2*con_type+1


class GaussianOrbBasis(object):
    def __init__(self, centers, shell_map, num_exponents, num_contractions, con_types, con_coeffs, exponents):
        """
           **Arguments:**

           centers
                A numpy array with shape (N, 3) with Cartesian centers.

           shell_map
                An array with the center index for each shell.

           num_exponents
                The number of exponents in each shell.

           num_contractions
                The number of contractions in each shell. This is used to
                implement optimized general contractions.

           con_types
                An array with contraction types: 0 = S, 1 = P, 2 = Cartesian D,
                3 = Cartesian F, ..., -2 = pure D, -3 = pure F, ...
                One contraction type is present for each contraction in each
                shell. The so-called SP type is implemented as a shell
                with two contractions, one of type S and one of type P.

           con_coeffs
                The contraction coefficients of the primitives for each
                contraction in each shell.

           exponents
                The exponents of the primitives in one shell.

           The number of primitives in shell i is num_exponents[i]*num_contractions[i].

           Convention for basis functions of a given contraction type:

           The order of the pure shells is based on the order of real spherical
           harmonics: http://en.wikipedia.org/wiki/Table_of_spherical_harmonics
           First the +- linear combination of highest angular momentum, then
           the ++ combination of highest angular momentum, keep repeating and
           finally take angular momention zero (without making a linear
           combination). The order of the Cartesian shells is sorted
           alhpabetically. The SP shell type is S first, then P. Some examples:

           shell_type=0, S:
             0 -> 1
           shell_type=1, P:
             0 -> x
             1 -> y
             2 -> z
           shell_type=2, Cartesian D:
             0 -> xx
             1 -> xy
             2 -> xz
             3 -> yy
             4 -> yz
             5 -> zz
           shell_type=3, Cartesian F:
             0 -> xxx
             1 -> xxy
             2 -> xxz
             3 -> xyy
             4 -> xyz
             5 -> xzz
             6 -> yyy
             7 -> yyz
             8 -> yzz
             9 -> zzz
           shell_type=-1, SP:
             0 -> 1
             1 -> x
             2 -> y
             3 -> z
           shell_type=-2, pure D:
             0 -> zz
             1 -> yz
             2 -> xz
             3 -> xx-yy
             4 -> xy
           shell_type=-3, pure F:
             6 -> zzz
             5 -> yzz
             4 -> xzz
             3 -> xxz-yyz
             2 -> xyz
             1 -> 3xxy-yyy
             0 -> xxx-3xyy
        """
        # This is the only thing that one may modify in place externally:
        self._centers = centers
        # All other fields are stored as internal parameters. Once they are set,
        # they are no supposed to be modified.
        self._shell_map = shell_map
        self._num_exponents = num_exponents
        self._num_contractions = num_contractions
        self._con_types = con_types
        self._con_coeffs = con_coeffs
        self._exponents = exponents
        # derived property, read only
        self._nbasis = sum(get_con_nbasis(con_type) for con_type in con_types)

    def get_centers(self):
        return self._centers.view()

    centers = property(get_centers)

    def get_nshell(self):
        return len(self._shell_map)

    nshell = property(get_nshell)

    def get_nbasis(self):
        return self._nbasis

    nbasis = property(get_nbasis)

    @classmethod
    def from_fchk(cls, fchk):
        shell_types = fchk.fields["Shell types"]
        shell_map = fchk.fields["Shell to atom map"] - 1
        num_exponents = fchk.fields["Number of primitives per shell"]
        ccoeffs_level1 = fchk.fields["Contraction coefficients"]
        ccoeffs_level2 = fchk.fields.get("P(S=P) Contraction coefficients")
        exponents = fchk.fields["Primitive exponents"]

        num_contractions = []
        con_types = []
        con_coeffs = []
        counter = 0
        for i, n in enumerate(num_exponents):
            if shell_types[i] == -1:
                # Special treatment for SP shell type
                num_contractions.append(2)
                con_types.append(0)
                con_types.append(1)
                tmp = np.array([
                    ccoeffs_level1[counter:counter+n],
                    ccoeffs_level2[counter:counter+n]
                ])
                con_coeffs.append(tmp.transpose().ravel())
            else:
                num_contractions.append(1)
                con_types.append(shell_types[i])
                con_coeffs.append(ccoeffs_level1[counter:counter+n])
            counter += n
        num_contractions = np.array(num_contractions)
        con_types = np.array(con_types)
        con_coeffs = np.concatenate(con_coeffs)


        result = cls(fchk.molecule.coordinates, shell_map, num_exponents,
                     num_contractions, con_types, con_coeffs, exponents)

        # permutation of the basis functions
        g03_reordering = {
          -6: np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
          -5: np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
          -4: np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]),
          -3: np.array([0, 1, 2, 3, 4, 5, 6]),
          -2: np.array([0, 1, 2, 3, 4]),
          -1: np.array([0, 1, 2, 3]),
           0: np.array([0]),
           1: np.array([0, 1, 2]),
           2: np.array([0, 3, 4, 1, 5, 2]),
           3: np.array([0, 4, 5, 3, 9, 6, 1, 8, 7, 2]),
           4: np.arange(15)[::-1],
           5: np.arange(21)[::-1],
           6: np.arange(28)[::-1],
        }
        offset = 0
        permutation = []
        for shell_type in shell_types:
            permutation.extend(g03_reordering[shell_type]+len(permutation))
        result.g03_permutation = np.array(permutation, dtype=int)

        return result
