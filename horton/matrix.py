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
"""Two- and four-dimensional matrix implementations

   The purpose of this module is to provide a generic API for different
   implementations of real-valued double precision matrix storage and
   operations.

   Two-dimensional matrices are supposed to be symmetric and are used to
   represent one-body operators and 1DRDMs. Four-dimensional matrices are used
   to represent two-body operators, which are invariant under the following
   interchanges of indexes::

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

from horton.log import log


__all__ = [
    'LinalgFactory', 'LinalgObject', 'Expansion', 'OneBody',
    'DenseLinalgFactory', 'DenseExpansion', 'DenseOneBody', 'DenseTwoBody',
]


class LinalgFactory(object):
    """A collection of compatible matrix and linear algebra routines.

       This is just an abstract base class that serves as a template for
       specific implementations.
    """
    def __init__(self, default_nbasis=None):
        '''
           **Optional arguments:**

           default_nbasis
                The default basis size when constructing new
                operators/expansions.
        '''
        self.default_nbasis = default_nbasis

    def create_expansion(self, nbasis=None):
        raise NotImplementedError

    def create_one_body(self, nbasis=None):
        raise NotImplementedError

    def create_two_body(self, nbasis=None):
        raise NotImplementedError

    def error_eigen(self, ham, overlap, expansion, epsilons):
        raise NotImplementedError

    def diagonalize(self, ham, overlap, expansion, epsilons):
        raise NotImplementedError

    def get_memory_one_body(self, nbasis=None):
        raise NotImplementedError

    def get_memory_two_body(self, nbasis=None):
        raise NotImplementedError


class LinalgObject(object):
    def apply_basis_permutation(self, permutation):
        raise NotImplementedError

    def apply_basis_signs(self, signs):
        raise NotImplementedError

    @classmethod
    def from_hdf5(cls, grp, lf):
        raise NotImplementedError

    def to_hdf5(self, grp):
        raise NotImplementedError

    def __clear__(self):
        self.clear()

    def clear(self):
        raise NotImplementedError

    def copy(self):
        raise NotImplementedError

    def assign(self, other):
        raise NotImplementedError


class Expansion(LinalgObject):
    def __init__(self, nbasis, nfn=None):
        raise NotImplementedError

    def check_normalization(self, olp, eps=1e-4):
        raise NotImplementedError


class OneBody(LinalgObject):
    def __init__(self, nbasis):
        raise NotImplementedError

    def set_element(self, i, j, value):
        raise NotImplementedError

    def get_element(self, i, j):
        raise NotImplementedError

    def iadd(self, other, factor=1):
        raise NotImplementedError

    def expectation_value(self, dm):
        raise NotImplementedError

    def trace(self):
        raise NotImplementedError

    def itranspose(self):
        raise NotImplementedError

    def iscale(self, factor):
        raise NotImplementedError

    def dot(self, vec0, vec1):
        raise NotImplementedError


class DenseLinalgFactory(LinalgFactory):
    def create_expansion(self, nbasis=None, nfn=None):
        nbasis = nbasis or self.default_nbasis
        return DenseExpansion(nbasis, nfn)

    def _check_expansion_init_args(self, expansion, nbasis=None, nfn=None):
        nbasis = nbasis or self.default_nbasis
        expansion.__check_init_args__(nbasis, nfn)

    create_expansion.__check_init_args__ = _check_expansion_init_args


    def create_one_body(self, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        return DenseOneBody(nbasis)

    def _check_one_body_init_args(self, one_body, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        one_body.__check_init_args__(nbasis)

    create_one_body.__check_init_args__ = _check_one_body_init_args


    def create_two_body(self, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        return DenseTwoBody(nbasis)

    def _check_two_body_init_args(self, two_body, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        two_body.__check_init_args__(nbasis)

    create_two_body.__check_init_args__ = _check_two_body_init_args


    @staticmethod
    def error_eigen(fock, overlap, expansion):
        """Compute the error of the orbitals with respect to the eigenproblem

           **Arguments:**

           fock
                A DenseOneBody Hamiltonian (or Fock) operator.

           overlap
                A DenseOneBody overlap operator.

           expansion
                An expansion object containing the current orbitals/eginvectors.
        """
        errors = np.dot(fock._array, expansion.coeffs) \
                 - expansion.energies*np.dot(overlap._array, expansion.coeffs)
        return np.sqrt((errors**2).mean())


    @staticmethod
    def diagonalize(fock, overlap=None):
        """Generalized eigen solver for the given Hamiltonian and overlap.

           **Arguments:**

           fock
                A DenseOneBody Hamiltonian (or Fock) operator.

           overlap
                A DenseOneBody overlap operator.

        """
        from scipy.linalg import eigh
        if overlap is None:
            return eigh(fock._array)
        else:
            return eigh(fock._array, overlap._array)

    def get_memory_one_body(self, nbasis=None):
        return nbasis**2*8

    def get_memory_two_body(self, nbasis=None):
        return nbasis**4*8


class DenseExpansion(Expansion):
    """An expansion of several functions in a basis with a dense matrix of
       coefficients. The implementation is such that the columns of self._array
       contain the orbitals.
    """
    def __init__(self, nbasis, nfn=None):
        """
           **Arguments:**

           nbasis
                The number of basis functions.

           **Optional arguments:**

           nfn
                The number of functions to store. Defaults to nbasis.

           do_energies
                Also allocate an array to store an energy corresponding to each
                function.
        """
        if nfn is None:
            nfn = nbasis
        self._coeffs = np.zeros((nbasis, nfn), float)
        self._energies = np.zeros(nfn, float)
        self._occupations = np.zeros(nfn, float)
        log.mem.announce(self._coeffs.nbytes + self._energies.nbytes + self._occupations.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_coeffs') and hasattr(self, '_energies') and hasattr(self, '_occupations'):
                log.mem.denounce(self._coeffs.nbytes + self._energies.nbytes + self._occupations.nbytes)

    def __check_init_args__(self, nbasis, nfn=None):
        if nfn is None:
            nfn = nbasis
        assert nbasis == self.nbasis
        assert nfn == self.nfn

    @classmethod
    def from_hdf5(cls, grp, lf):
        if grp.attrs['class'] != cls.__name__:
            raise TypeError('The class of the expansion in the HDF5 file does not match.')
        nbasis, nfn = grp['coeffs'].shape
        result = cls(nbasis, nfn)
        grp['coeffs'].read_direct(result._coeffs)
        grp['energies'].read_direct(result._energies)
        grp['occupations'].read_direct(result._occupations)
        return result

    def to_hdf5(self, grp):
        grp.attrs['class'] = self.__class__.__name__
        grp['coeffs'] = self._coeffs
        grp['energies'] = self._energies
        grp['occupations'] = self._occupations

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._coeffs.shape[0]

    nbasis = property(_get_nbasis)

    def _get_nfn(self):
        '''The number of orbitals (or functions in general)'''
        return self._coeffs.shape[1]

    nfn = property(_get_nfn)

    def _get_coeffs(self):
        '''The matrix with the expansion coefficients'''
        return self._coeffs.view()

    coeffs = property(_get_coeffs)

    def _get_energies(self):
        '''The orbital energies'''
        return self._energies

    energies = property(_get_energies)

    def _get_occupations(self):
        '''The orbital occupations'''
        return self._occupations

    occupations = property(_get_occupations)

    def clear(self):
        self._coeffs[:] = 0.0
        self._energies[:] = 0.0
        self._occupations[:] = 0.0

    def copy(self):
        result = DenseExpansion(self.nbasis, self.nfn)
        result._coeffs[:] = self._coeffs
        result._energies[:] = self._energies
        result._occupations[:] = self._occupations
        return result

    def check_normalization(self, olp, eps=1e-4):
        '''Run an internal test to see if the orbitals are normalized

           **Arguments:**

           olp
                The overlap one_body operators

           **Optional arguments:**

           eps
                The allowed deviation from unity, very loose by default.
        '''
        for i in xrange(self.nfn):
            if self.occupations[i] == 0:
                continue
            norm = olp.dot(self._coeffs[:,i], self._coeffs[:,i])
            #print i, norm
            assert abs(norm-1) < eps, 'The orbitals are not normalized!'

    def to_dm(self, output=None, factor=None):
        """Compute the density matrix

           **Optional arguments:**

           output
                An output density matrix (DenseOneBody instance).

           factor
                When given, the density matrix is added with the given prefactor
                to the output argument. If not given, the original contents of
                dm are overwritten. This argument implies that the dm output
                argument must also be present.
        """
        # parse first argument
        if output is None:
            dm = DenseOneBody(self.nbasis)
            if factor is not None:
                raise TypeError('When the factor argument is given, the output argument must be a density matrix.')
        else:
            dm = output
        if factor is None:
            dm._array[:] = np.dot(self._coeffs*self.occupations, self._coeffs.T)
        else:
            dm._array[:] += factor*np.dot(self._coeffs*self.occupations, self._coeffs.T)
        return dm

    def from_fock(self, fock, overlap):
        '''Diagonalize a Fock matrix to obtain orbitals and energies'''
        evals, evecs = DenseLinalgFactory.diagonalize(fock, overlap)
        self._energies[:] = evals[:self.nfn]
        self._coeffs[:] = evecs[:,:self.nfn]

    def derive_naturals(self, dm, overlap):
        '''
           **Arguments**:

           dm
                A DenseOneBody object with the density matrix

           overlap
                A DenseOneBody object with the overlap matrix

           **Optional arguments:**
        '''
        # Construct a level-shifted operator
        occ = overlap.copy()
        occ.idot(dm)
        occ.idot(overlap)
        # diagonalize and compute eigenvalues
        evals, evecs = DenseLinalgFactory.diagonalize(occ, overlap)
        self._coeffs[:] = evecs[:,:self.nfn]
        self._occupations[:] = evals
        self._energies[:] = 0.0

    def apply_basis_permutation(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        self._coeffs[:] = self.coeffs[permutation]

    def apply_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.'''
        self._coeffs *= signs.reshape(-1,1)

    def assign(self, other):
        if not isinstance(other, DenseExpansion):
            raise TypeError('The other object must also be DenseExpansion instance.')
        self._coeffs[:] = other._coeffs
        self._energies[:] = other._energies
        self._occupations[:] = other._occupations

    def get_homo_index(self, offset=0):
        '''Return the index of a HOMO orbital.'''
        if offset < 0:
            raise ValueError('Offset must be zero or positive.')
        homo_indexes = self.occupations.nonzero()[0]
        if len(homo_indexes) > offset:
            return homo_indexes[len(homo_indexes)-offset-1]

    def get_homo_energy(self, offset=0):
        '''Return a homo energy

           **Optional arguments**:

           offset
                By default, the (highest) homo energy is returned. When this
                index is above zero, the corresponding lower homo energy is
                returned.
        '''
        index = self.get_homo_index(offset)
        if index is not None:
            return self.energies[index]

    homo_energy = property(get_homo_energy)

    def get_lumo_index(self, offset=0):
        '''Return the index of a LUMO orbital.'''
        if offset < 0:
            raise ValueError('Offset must be zero or positive.')
        lumo_indexes = (self.occupations==0.0).nonzero()[0]
        if len(lumo_indexes) > offset:
            return lumo_indexes[offset]

    def get_lumo_energy(self, offset=0):
        '''Return a lumo energy

           **Optional arguments**:

           offset
                By default, the (lowest) lumo energy is returned. When this
                index is above zero, the corresponding higher homo energy is
                returned.
        '''
        index = self.get_lumo_index(offset)
        if index is not None:
            return self.energies[index]

    lumo_energy = property(get_lumo_energy)


class DenseOneBody(OneBody):
    """Dense symmetric two-dimensional matrix, also used for density matrices.

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
        self._array = np.zeros((nbasis, nbasis), float)
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)

    def __check_init_args__(self, nbasis):
        assert nbasis == self.nbasis

    @classmethod
    def from_hdf5(cls, grp, lf):
        nbasis = grp['array'].shape[0]
        result = cls(nbasis)
        grp['array'].read_direct(result._array)
        return result

    def to_hdf5(self, grp):
        grp.attrs['class'] = self.__class__.__name__
        grp['array'] = self._array

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._array.shape[0]

    nbasis = property(_get_nbasis)

    def set_element(self, i, j, value):
        self._array[i,j] = value
        self._array[j,i] = value

    def get_element(self, i, j):
        return self._array[i,j]

    def assign(self, other):
        if not isinstance(other, DenseOneBody):
            raise TypeError('The other object must also be DenseOneBody instance.')
        self._array[:] = other._array

    def set_diagonal(self, value):
        '''Set diagonal elements to value'''
        np.fill_diagonal(self._array, value)

    def copy(self):
        '''Return a copy of the current one-body operator'''
        result = DenseOneBody(self.nbasis)
        result._array[:] = self._array
        return result

    def new(self):
        '''Return a new one-body operator with the same nbasis'''
        return DenseOneBody(self.nbasis)

    def _check_new_init_args(self, one_body, nbasis=None):
        one_body.__check_init_args__(self.nbasis)

    new.__check_init_args__ = _check_new_init_args


    def check_symmetry(self):
        '''Check the symmetry of the array. For testing only.'''
        assert abs(self._array - self._array.T).max() == 0.0

    def clear(self):
        '''Resets array to zeros element-wise.'''
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

    def iscale(self, factor):
        self._array *= factor

    def dot(self, vec0, vec1):
        return np.dot(vec0, np.dot(self._array, vec1))

    def idot(self, other):
        '''Inplace dot operator'''
        self._array[:] = np.dot(self._array, other._array)

    def distance(self, other):
        '''Maximum difference between self and other one body object'''
        return abs(self._array.ravel() - other._array.ravel()).max()

    def apply_basis_permutation(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        self._array[:] = self._array[permutation]
        self._array[:] = self._array[:,permutation]

    def apply_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.'''
        self._array *= signs
        self._array *= signs.reshape(-1,1)


class DenseTwoBody(LinalgObject):
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
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)

    def __check_init_args__(self, nbasis):
        assert nbasis == self.nbasis

    @classmethod
    def from_hdf5(cls, grp, lf):
        nbasis = grp['array'].shape[0]
        result = cls(nbasis)
        grp['array'].read_direct(result._array)
        return result

    def to_hdf5(self, grp):
        grp.attrs['class'] = self.__class__.__name__
        grp['array'] = self._array

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._array.shape[0]

    nbasis = property(_get_nbasis)

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
        output._array[:] = np.tensordot(self._array, dm._array, ([1,3], [1,0]))

    def apply_exchange(self, dm, output):
        """Compute the exchange dot product with a density matrix."""
        if not isinstance(dm, DenseOneBody):
            raise TypeError('The dm argument must be a DenseOneBody class')
        if not isinstance(output, DenseOneBody):
            raise TypeError('The output argument must be a DenseOneBody class')
        output._array[:] = np.tensordot(self._array, dm._array, ([1,2], [1,0]))

    def clear(self):
        self._array[:] = 0.0

    def apply_basis_permutation(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        self._array[:] = self._array[permutation]
        self._array[:] = self._array[:,permutation]
        self._array[:] = self._array[:,:,permutation]
        self._array[:] = self._array[:,:,:,permutation]

    def apply_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.'''
        self._array *= signs
        self._array *= signs.reshape(-1,1)
        self._array *= signs.reshape(-1,-1,1)
        self._array *= signs.reshape(-1,-1,-1,1)
