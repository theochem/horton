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
"""Two- and four-dimensional matrix implementations.

   The purpose of this module is to provide a generic API for different
   implementations of real-valued double precision matrix storage and
   operations.

   Two-dimensional matrices are supposed to be symmetric and are used to
   represent one-body operators and 1DRDMs. Four-dimensional matrices are used
   to represent two-body operators, which are invariant under the following
   interchanges of indexes:
            <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
            <il|kj> = <jk|li> = <kj|il> = <li|jk>
   This module assumes physicists notation for the two-particle operators. It is
   up to the specific implementations of the matrices to make use of these
   symmetries.

   One should use these matrix implementations without accessing the internals
   of each class, i.e. without accessing attributes or methodsthat start with an
   underscore.

   In order to avoid temporaries when working with these arrays, all operations
   are defined as in-place operations. This forces the user to allocate all
   memory in advance, which can then be moved out of the inner loops.
"""


import numpy as np


__all__ = [
    'DenseExpansion', 'DenseOneBody', 'DenseTwoBody',
    'error_eigen', 'diagonalize'
]



class DenseExpansion(object):
    """An expansion of several functions in a basis with a dense matrix of
       coefficients.
    """
    def __init__(self, nfn, nbasis):
        self._coeffs = np.zeros((nbasis, nfn), float)

    def get_nbasis(self):
        return self._coeffs.shape[1]

    nbasis = property(get_nbasis)

    def get_density_matrix(self, noc):
        """Get the density matrix

           **Arguments:**

           noc
                The number of 'occupied' functions
        """
        result = DenseOneBody(self.nbasis)
        occupied = self._coeffs[:,:noc]
        result._array[:] = np.dot(occupied, occupied.T)
        return result


class DenseOneBody(object):
    """Dense symmetric two-dimensional matrix.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """
    def __init__(self, size=None):
        """
           **Arguments:**

           size
                The number of basis functions.
        """
        self._array = np.zeros((size, size), float)

    def get_size(self):
        return self._array.shape[0]

    size = property(get_size)

    def set_element(self, i, j, value):
        self._array[i,j] = value
        self._array[j,i] = value

    def get_element(self, i, j):
        return self._array[i,j]

    def check_symmetry(self):
        """Check the symmetry of the array."""
        assert abs(self._array - self._array.T).max() == 0.0

    def reset(self):
        self._array[:] = 0.0

    def iadd(self, other, factor=1):
        self._array += other._array*factor

    def expectation_value(self, dm):
        return np.dot(self._array.ravel(), dm._array.ravel())


class DenseTwoBody(object):
    """Dense symmetric four-dimensional matrix.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """
    def __init__(self, size):
        """
           **Arguments:**

           size
                The number of basis functions.
        """
        self._array = np.zeros((size, size, size, size), float)

    def get_size(self):
        return self._array.shape[0]

    size = property(get_size)

    def set_element(self, i, j, k, l, value):
        #    <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
        #    <il|kj> = <jk|li> = <kj|il> = <li|jk>
        self._array[i,j,k,l] = value
        self._array[j,i,l,k] = value
        self._array[k,l,i,j] = value
        self._array[l,k,j,i] = value
        self._array[i,l,k,j] = value
        self._array[j,k,l,i] = value
        self._array[k,j,i,l] = value
        self._array[l,i,j,k] = value

    def get_element(self, i, j, k, l):
        return self._array[i,j, k, l]

    def check_symmetry(self):
        """Check the symmetry of the array."""
        assert abs(self._array - self._array.transpose(1,0,3,2)).max() == 0.0
        assert abs(self._array - self._array.transpose(2,3,0,1)).max() == 0.0
        assert abs(self._array - self._array.transpose(3,2,1,0)).max() == 0.0
        assert abs(self._array - self._array.transpose(2,1,0,3)).max() == 0.0
        assert abs(self._array - self._array.transpose(3,0,1,2)).max() == 0.0
        assert abs(self._array - self._array.transpose(0,3,2,1)).max() == 0.0
        assert abs(self._array - self._array.transpose(1,2,3,0)).max() == 0.0

    def apply_direct(self, dm, output):
        """Compute the direct dot product with a density matrix."""
        if not isinstance(dm, DenseOneBody):
            raise TypeError('The dm argument must be a DenseOneBody class')
        if not isinstance(output, DenseOneBody):
            raise TypeError('The output argument must be a DenseOneBody class')
        output._array[:] = np.tensordot(self._array, dm._array, ([1,3], [0,1]))

    def apply_exchange(self, dm, output):
        """Compute the exchange dot product with a density matrix."""
        if not isinstance(dm, DenseOneBody):
            raise TypeError('The dm argument must be a DenseOneBody class')
        if not isinstance(output, DenseOneBody):
            raise TypeError('The output argument must be a DenseOneBody class')
        output._array[:] = np.tensordot(self._array, dm._array, ([1,2], [0,1]))


def error_eigen(ham, overlap, wfn):
    errors = np.dot(ham._array, wfn._expansion._coeffs) \
             - wfn._epsilons*np.dot(overlap._array, wfn._expansion._coeffs)
    return np.sqrt((errors**2).mean())


def diagonalize(ham, overlap, wfn):
    from scipy.linalg import eigh
    evals, evecs = eigh(ham._array, overlap._array)
    wfn._expansion._coeffs[:] = evecs
    wfn._epsilons[:] = evals
