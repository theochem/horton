# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
   of each class, i.e. without accessing attributes or methods that start with
   an underscore.

   In order to avoid temporaries when working with arrays, the methods do
   not return arrays. Instead such methods are an in place operation or have
   output arguments. This forces the user to allocate all memory in advance,
   which can then be moved out of the loops. The initial implementation (the
   Dense... classes) are just a proof of concept and may therefore contain
   internals that still make temporaries. This fixed later with an alternative
   implementation.
"""


import numpy as np


__all__ = [
    'LinalgFactory', 'DenseLinalgFactory', 'DenseExpansion', 'DenseOneBody',
    'DenseTwoBody',
]


class LinalgFactory(object):
    """A collection of compatible matrix and linear algebra routines.

       This is just an abstract base class that serves as a template for
       specific implementations.
    """
    def create_expansion(self, nbasis):
        raise NotImplementedError

    def create_one_body(self, nbasis):
        raise NotImplementedError

    def create_two_body(self, nbasis):
        raise NotImplementedError

    def error_eigen(self, ham, overlap, expansion, epsilons):
        raise NotImplementedError

    def diagonalize(self, ham, overlap, expansion, epsilons):
        raise NotImplementedError

    def get_memory_one_body(self, nbasis):
        raise NotImplementedError

    def get_memory_two_body(self, nbasis):
        raise NotImplementedError


class DenseLinalgFactory(LinalgFactory):
    def create_expansion(self, nbasis, nfn, do_energies=False):
        return DenseExpansion(nbasis, nfn, do_energies)

    def create_one_body(self, nbasis):
        return DenseOneBody(nbasis)

    def create_two_body(self, nbasis):
        return DenseTwoBody(nbasis)

    def error_eigen(self, ham, overlap, expansion):
        """Compute the error of the orbitals with respect to the eigenproblem

           **Arguments:**

           ham
                A DenseOneBody Hamiltonian (or Fock) operator.

           overlap
                A DenseOneBody overlap operator.

           expansion
                An expansion object containing the current orbitals/eginvectors.

           epsilons
                An array with the orbital energies.
        """
        errors = np.dot(ham._array, expansion.coeffs) \
                 - expansion.energies*np.dot(overlap._array, expansion.coeffs)
        return np.sqrt((errors**2).mean())


    def diagonalize(self, ham, overlap, expansion):
        """Generalized eigen solver for the given Hamiltonian and overlap.

           **Arguments:**

           ham
                A DenseOneBody Hamiltonian (or Fock) operator.

           overlap
                A DenseOneBody overlap operator.

        """
        from scipy.linalg import eigh
        evals, evecs = eigh(ham._array, overlap._array)
        expansion.coeffs[:] = evecs
        expansion.energies[:] = evals

    def get_memory_one_body(self, nbasis):
        return nbasis**2*8

    def get_memory_two_body(self, nbasis):
        return nbasis**4*8


class DenseExpansion(object):
    """An expansion of several functions in a basis with a dense matrix of
       coefficients. The implementation is such that the columns of self._array
       contain the orbitals
    """
    def __init__(self, nbasis, nfn, do_energies=False):
        """
           **Arguments:**

           nbasis
                The number of basis functions

           nfn
                The number of functions to store

           **Optional arguments:**

           do_energies
                Also allocate an array to store an energy corresponding to each
                function.
        """
        self._coeffs = np.zeros((nbasis, nfn), float)
        if do_energies:
            self._energies = np.zeros(nfn, float)
        else:
            self._energies = None

    def get_nbasis(self):
        return self._coeffs.shape[1]

    nbasis = property(get_nbasis)

    def get_coeffs(self):
        return self._coeffs.view()

    coeffs = property(get_coeffs)

    def get_energies(self):
        return self._energies.view()

    energies = property(get_energies)

    def compute_density_matrix(self, noc, dm, factor=None):
        """Compute the density matrix

           **Arguments:**

           noc
                The number of 'occupied' functions

           dm
                An output density matrix. This must be a DenseOneBody instance.

           **Optional arguments:**

           factor
                When given, the density matrix is added with the given prefactor
                to the output argument. If not given, the original contents of
                dm are overwritten.
        """
        result = DenseOneBody(self.nbasis)
        occupied = self._coeffs[:,:noc]
        if factor is None:
            dm._array[:] = np.dot(occupied, occupied.T)
        else:
            dm._array[:] += factor*np.dot(occupied, occupied.T)

    def apply_basis_permutation(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        self._coeffs[:] = self.coeffs[permutation]


class DenseOneBody(object):
    """Dense symmetric two-dimensional matrix, also used for density matrices.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """
    def __init__(self, nbasis=None):
        """
           **Arguments:**

           nbasis
                The number of basis functions.
        """
        self._array = np.zeros((nbasis, nbasis), float)

    def get_nbasis(self):
        return self._array.shape[0]

    nbasis = property(get_nbasis)

    def set_element(self, i, j, value):
        self._array[i,j] = value
        self._array[j,i] = value

    def get_element(self, i, j):
        return self._array[i,j]

    def check_symmetry(self):
        '''Check the symmetry of the array. For testing only.'''
        assert abs(self._array - self._array.T).max() == 0.0

    def reset(self):
        self._array[:] = 0.0

    def iadd(self, other, factor=1):
        self._array += other._array*factor

    def expectation_value(self, dm):
        return np.dot(self._array.ravel(), dm._array.ravel())

    def trace(self):
        return np.trace(self._array)

    def itranspose(self):
        '''In-place transpose'''
        self._array = self._array.T

    def dot(self, vec0, vec1):
        return np.dot(vec0, np.dot(self._array, vec1))

    def apply_basis_permutation(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        self._array[:] = self._array[permutation]
        self._array[:] = self._array[:,permutation]


class DenseTwoBody(object):
    """Dense symmetric four-dimensional matrix.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """
    def __init__(self, nbasis):
        """
           **Arguments:**

           nbasis
                The number of basis functions.
        """
        self._array = np.zeros((nbasis, nbasis, nbasis, nbasis), float)

    def get_nbasis(self):
        return self._array.shape[0]

    nbasis = property(get_nbasis)

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

    def apply_basis_permutation(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        self._array[:] = self._array[permutation]
        self._array[:] = self._array[:,permutation]
        self._array[:] = self._array[:,:,permutation]
        self._array[:] = self._array[:,:,:,permutation]
