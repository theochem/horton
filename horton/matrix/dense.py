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
"""Dense matrix implementations"""


import numpy as np
from horton.log import log
from horton.matrix.base import LinalgFactory, LinalgObject, Expansion, OneIndex, TwoIndex, ThreeIndex, FourIndex


__all__ = [
    'DenseLinalgFactory', 'DenseExpansion',
    'DenseOneIndex', 'DenseTwoIndex', 'DenseThreeIndex', 'DenseFourIndex'
]


class DenseLinalgFactory(LinalgFactory):
    def create_expansion(self, nbasis=None, nfn=None):
        nbasis = nbasis or self.default_nbasis
        nfn = nfn or self.default_nbasis
        return DenseExpansion(nbasis, nfn)

    def _check_expansion_init_args(self, expansion, nbasis=None, nfn=None):
        nbasis = nbasis or self.default_nbasis
        nfn = nfn or self.default_nbasis
        expansion.__check_init_args__(nbasis, nfn)

    create_expansion.__check_init_args__ = _check_expansion_init_args

    def create_two_index(self, nbasis=None, nfn=None):
        nbasis = nbasis or self.default_nbasis
        nfn = nfn or self.default_nbasis
        return DenseTwoIndex(nbasis, nfn)

    def _check_two_index_init_args(self, two_index, nbasis=None, nfn=None):
        nbasis = nbasis or self.default_nbasis
        nfn = nfn or self.default_nbasis
        two_index.__check_init_args__(nbasis)

    create_two_index.__check_init_args__ = _check_two_index_init_args

    def create_one_index(self, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        return DenseOneIndex(nbasis)

    def _check_one_index_init_args(self, two_index, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        two_index.__check_init_args__(nbasis)

    create_one_index.__check_init_args__ = _check_one_index_init_args

    def create_four_index(self, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        return DenseFourIndex(nbasis)

    def _check_four_index_init_args(self, four_index, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        four_index.__check_init_args__(nbasis)

    create_four_index.__check_init_args__ = _check_four_index_init_args

    def create_three_index(self, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        return DenseThreeIndex(nbasis)

    def _check_three_index_init_args(self, three_index, nbasis=None):
        nbasis = nbasis or self.default_nbasis
        three_index.__check_init_args__(nbasis)

    create_three_index.__check_init_args__ = _check_three_index_init_args

    @staticmethod
    def error_eigen(fock, overlap, expansion):
        """Compute the error of the orbitals with respect to the eigenproblem

           **Arguments:**

           fock
                A DenseTwoIndex Hamiltonian (or Fock) operator.

           overlap
                A DenseTwoIndex overlap operator.

           expansion
                An expansion object containing the current orbitals/eginvectors.
        """
        errors = np.dot(fock._get_dense(), (expansion.coeffs)) \
                 - expansion.energies*np.dot(overlap._get_dense(), (expansion.coeffs))
        return np.sqrt((abs(errors)**2).mean())


    @staticmethod
    def diagonalize(fock, overlap=None):
        """Generalized eigen solver for the given Hamiltonian and overlap.

           **Arguments:**

           fock
                A DenseTwoIndex Hamiltonian (or Fock) operator.

           overlap
                A DenseTwoIndex overlap operator.

        """
        from scipy.linalg import eigh
        if overlap is None:
            return eigh(fock._get_dense())
        else:
            return eigh(fock._get_dense(), overlap._get_dense())

    def get_memory_two_index(self, nbasis=None):
        return nbasis**2*8

    def get_memory_four_index(self, nbasis=None):
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
        self._coeffs = np.zeros((nbasis, nfn))
        self._energies = np.zeros(nfn)
        self._occupations = np.zeros(nfn)
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
    def from_hdf5(cls, grp):
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

    def copy_occ_model(self, other):
        self._occupations[:] = other._occupations[:]

    def dotarray(self, other, array):
        '''Dot product of other DenseExpansion and numpy array'''
        self._coeffs[:] = np.dot(other.coeffs, array)

    def check_normalization(self, olp, eps=1e-4):
        '''Run an internal test to see if the orbitals are normalized

           **Arguments:**

           olp
                The overlap two_index operators

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


    def to_dm(self, output=None, factor=None, other=None):
        """Compute the density matrix

           **Optional arguments:**

           output
                An output density matrix (DenseTwoIndex instance).

           factor
                When given, the density matrix is added with the given prefactor
                to the output argument. If not given, the original contents of
                dm are overwritten. This argument implies that the dm output
                argument must also be present.
        """
        # parse first argument
        if output is None:
            dm = DenseTwoIndex(self.nbasis, self.nbasis)
            if factor is not None:
                raise TypeError('When the factor argument is given, the output argument must be a density matrix.')
        else:
            dm = output
        if factor is None:
            if other is None:
                dm.assign_array(np.dot(self._coeffs*self.occupations, self._coeffs.T))
            else:
                dm.assign_array(np.dot(self._coeffs*self.occupations, other._coeffs.T))
        else:
            if other is None:
                dm._iadd_array(factor*np.dot(self._coeffs*self.occupations, self._coeffs.T))
            else:
                dm._iadd_array(factor*np.dot(self._coeffs*self.occupations, other._coeffs.T))
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
                A DenseTwoIndex object with the density matrix

           overlap
                A DenseTwoIndex object with the overlap matrix

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
            raise TypeError('The other object must also be DenseExpansion instance. Got ', type(other), ' instead.')
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

    def unitary_rotation(self, n=None):
        '''Apply random unitary transformation distributed with Haar measure'''
        from scipy import linalg as linalg
        from scipy import randn, diagonal, absolute, multiply
        if n is None:
            n = self.nbasis
        z = randn(n, n)
        q, r = linalg.qr(z)
        d = diagonal(r)
        ph = d/absolute(d)
        q = multiply(q, ph, q)
        self.coeffs[:] = np.dot(self.coeffs, q)

    def mixing(self, angle=45, index1=None, index2=None):
        '''Mix Homo and Lumo orbitals with angle 'angle'.
           This works only for singlets.
        '''
        if index2 == None:
            index2 = self.get_lumo_index()
        if index1 == None:
            index1 = self.get_homo_index()
        print index1, index2
        mixing = self.compute_rotation(index1, index2, angle)
        self.coeffs[:] = np.dot(self.coeffs, mixing)

    def compute_rotation(self, index1, index2, angle=45):
        '''Get Rotations matrix performing a Givens rotation.
        '''
        from scipy import math
        output = np.identity((self.nbasis))
        radians = (float(angle)/180.0)*math.pi
        output[index1,index1] = math.cos(radians)
        output[index2,index2] = math.cos(radians)
        if index1 < index2:
            output[index1,index2] = -math.sin(radians)
            output[index2,index1] = math.sin(radians)
        else:
            output[index1,index2] = math.sin(radians)
            output[index2,index1] = -math.sin(radians)
        return output

    def apply_swapping(self, swapvector, factor=1.0):
        '''Change order of MOs using pair-exchange'''
        for swapi in range(0,len(swapvector),2):
            index1 = swapvector[(swapi+1)]-1
            index2 = swapvector[swapi]-1
            if log.do_medium:
                log('  Swapping orbitals %i and %i' %(index1+1, index2+1))
            swapelement = np.array([index1, index2])
            neworder = np.array([index2, index1])
            self.coeffs[:,neworder] = self.coeffs[:,swapelement]*factor
            self.energies[neworder] = self.energies[swapelement]

    def apply_givens_rotation(self, givensrotation):
        '''Do Givens rotation
        '''
        for roti in range(0,len(givensrotation),3):
            index1 = givensrotation[roti]-1
            index2 = givensrotation[(roti+1)]-1
            angle = givensrotation[(roti+2)]
            if log.do_medium:
                log('  Rotating orbitals %i and %i with angle %3.1f'
                       %(index1+1, index2+1, angle))
            givens = self.compute_rotation(index1, index2, angle)
            self.coeffs[:] = np.dot(self.coeffs, givens)


class DenseOneIndex(OneIndex):
    """Dense one-dimensional matrix (vector), also used for (diagonal) density matrices.
    """
    def __init__(self, nbasis):
        """
           **Arguments:**

           nbasis
                The number of basis functions.
        """
        self._array = np.zeros((nbasis,1))
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            log.mem.denounce(self._array.nbytes)

    def __check_init_args__(self, nbasis):
        assert nbasis == self.nbasis

    @classmethod
    def from_hdf5(cls, grp):
        nbasis = grp['array'].shape[0]
        result = cls(nbasis)
        grp['array'].read_direct(result._array)
        return result

    def read_from_hdf5(self, grp):
        if grp.attrs['class'] != self.__class__.__name__:
            raise TypeError('The class of the one-index object in the HDF5 file does not match.')
        grp['array'].read_direct(self._array)

    def to_hdf5(self, grp):
        grp.attrs['class'] = self.__class__.__name__
        grp['array'] = self._array

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._array.shape[0]

    nbasis = property(_get_nbasis)

    def set_element(self, i, value):
        self._array[i] = value

    def get_element(self, i):
        return self._array[i]

    def set_view(self, ind1, ind2, value):
        self._array[ind1:ind2] = value

    def assign(self, other):
        if not isinstance(other, DenseOneIndex):
            raise TypeError('The other object must also be DenseTwoIndex instance. Got ', type(other), ' instead.')
        self._array[:] = other._array

    def copy(self):
        result = DenseOneIndex(self.nbasis)
        result._array[:] = self._array
        return result

    def copyview(self, ind1, ind2):
        dim = ind2-ind1
        result = DenseOneIndex(dim)
        result._array[:] = self._array[ind1:ind2]
        return result

    def clear(self):
        self._array[:] = 0.0

    def iadd(self, other, factor=1):
        self._array += other._array*factor

    def trace(self):
        return np.sum(self._array)

    def iscale(self, factor):
        self._array *= factor

    def dot(self, vec):
        return np.dot(self._array, vec)

    def idot(self, other):
        self._array[:] = self._array*other._array

    def distance(self, other):
        return np.linalg.norm((self._array.ravel() - other._array.ravel()))

    def contract_twoindex(self, other, select, factor=1.0, *args):
        '''Contract TwoIndex to OneIndex.

           **Arguments:**

           other
                A DenseTwoIndex object.

           select:
                '1': contract first index
                '2': contract second index
        '''
        ind = []
        for arg in args:
            ind.append(arg)
        if not ind:
            ind = [0, other._get_dense().shape[0], 0, other._get_dense().shape[1]]
        if select == 1:
            self._array[:] = (np.einsum('ab->b', other._get_dense()[ind[0]:ind[1], ind[2]:ind[3]])*factor)[np.newaxis].T
        elif select == 2:
            self._array[:] = (np.einsum('ab->a', other._get_dense()[ind[0]:ind[1], ind[2]:ind[3]])*factor)[np.newaxis].T
        else:
            raise NotImplementedError

    def contract_2twoindex(self, other, other2, select, factor=1.0):
        '''Contract two TwoIndex to OneIndex.

           **Arguments:**

           other/other2
                A DenseTwoIndex object.

           select:
                '1': contract first index
                '2': contract second index
        '''
        if select == 1:
            self._array[:] = (np.einsum('ab,ab->b', other._get_dense(), other2._get_dense())*factor)[np.newaxis].T
        elif select == 2:
            self._array[:] = (np.einsum('ab,ab->a', other._get_dense(), other2._get_dense())*factor)[np.newaxis].T
        else:
            raise NotImplementedError

    def assign_diagonal(self, other, factor=1.0, ind1=0, ind2=None):
        if not ind2:
            ind2 = other._get_dense().shape[0]
        self._array[:] = ((other._get_dense().diagonal())[ind1:ind2]*factor)[np.newaxis].T

    def compute_ps2_one_dm_ap1rog(self, geminal, lagrange, factor=1.0):
        '''Compute PS2 1-RDM for AP1roG
           This 1-DM is only used for the variational orbital optimization.
        '''
        nd_lagrange = lagrange._get_dense()
        nd_geminal = geminal._get_dense()

        nocc = nd_geminal.shape[0]
        value = 1+geminal.contract_twoindex(lagrange)

        self._array[:nocc] = value-np.einsum("ia,ia->i", nd_lagrange, nd_geminal)[np.newaxis].T
        self._array[nocc:] = np.einsum("ia,ia->a", nd_lagrange, nd_geminal)[np.newaxis].T
        self.iscale(factor)

    def compute_response_one_dm_ap1rog(self, geminal, lagrange, factor=1.0):
        '''Compute response 1-RDM for AP1roG'''
        nd_lagrange = lagrange._get_dense()
        nd_geminal = geminal._get_dense()

        nocc = nd_geminal.shape[0]
        self._array[:nocc] = 1-np.einsum("ia,ia->i", nd_lagrange, nd_geminal)[np.newaxis].T
        self._array[nocc:] = np.einsum("ia,ia->a", nd_lagrange, nd_geminal)[np.newaxis].T
        self.iscale(factor)


class DenseTwoIndex(TwoIndex):
    """Dense symmetric two-dimensional matrix, also used for density matrices.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """
    def __init__(self, nbasis, nfn=None):
        """
           **Arguments:**

           nbasis
                The number of basis functions.
        """
        if nfn is None:
            nfn = nbasis
        self._array = np.zeros((nbasis, nfn))
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)

    def __check_init_args__(self, nbasis):
        assert nbasis == self.nbasis

    @classmethod
    def from_hdf5(cls, grp):
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

    def _get_nfn(self):
        '''The number of orbitals (or functions in general)'''
        return self._array.shape[1]

    nfn = property(_get_nfn)

    def set_element(self, i, j, value):
        self._array[i,j] = value
        self._array[j,i] = value

    def get_element(self, i, j):
        return self._array[i,j]

    def _get_dense(self):
        '''Returns an immutable view of the current matrix.
            Not for use outside of the Matrix class.
        '''
        b = self._array.view()
        b.flags['WRITEABLE']=False
        return b

    def assign(self, other):
        if not isinstance(other, DenseTwoIndex):
            raise TypeError('The other object must also be DenseTwoIndex instance. Got ', type(other), ' instead.')
        self._array[:] = other._array

    def set_diagonal(self, value):
        '''Set diagonal elements to value'''
        if isinstance(value, DenseOneIndex):
            np.fill_diagonal(self._array, value._array.T)
        else:
            np.fill_diagonal(self._array, value)

    def assign_array(self, other, dim1=None, dim2=None):
        '''Assign underlying representation from a numpy array.

            ** Arguments: **
                other
                    A numpy ndarray
                dim1
                    Dimension 1 to reshape 'other'
                dim2
                    Dimension 2 to reshape 'other'
        '''
        if dim1:
            self._array[:] = other.reshape((dim1, dim2))
        else:
            self._array[:] = other

    def set_value(self, val):
        self._array[:] = val

    def _set_from_dense(self, arr):
        '''Set the underlying representation from a dense 2-dim array

            ** Arguments: **
                arr
                    A 2-dim numpy ndarray. This will be copied into the underlying representation.
        '''

        assert isinstance(arr, np.ndarray)
        self._array[:] = arr

    def copy(self):
        '''Return a copy of the current two-index operator'''
        result = DenseTwoIndex(self.nbasis, self.nfn)
        result._array[:] = self._array
        return result

    def copyview(self, ind1, ind2, ind3, ind4, factor=1.0):
        dim1 = ind2-ind1
        dim2 = ind4-ind3
        result = DenseTwoIndex(dim1, dim2)
        result._array[:] = self._array[ind1:ind2, ind3:ind4]*factor
        return result

    def assign_tril_fourindex(self, other, dim):
        tril = np.tril_indices(dim, -1)
        self._array[:] = (other._array[:,:,tril[0],tril[1]])[tril]

    def new(self):
        '''Return a new two-index operator with the same nbasis'''
        return DenseTwoIndex(self.nbasis, self.nbasis)

    def _check_new_init_args(self, two_index, nbasis=None):
        two_index.__check_init_args__(self.nbasis)

    new.__check_init_args__ = _check_new_init_args


    def check_symmetry(self):
        '''Check the symmetry of the array. For testing only.'''
        assert abs(self._array - self._array.T).max() == 0.0

    def isymmetrize(self):
        '''Symmetrize DenseTwoIndex using M_sym=(M+M^\dagger)/2'''
        self._array = (self._array+self._array.T)/2.0

    def symmetrize(self, factor=1.0):
        '''Symmetrize DenseTwoIndex using M_sym=(M+M^\dagger)/2'''
        result = DenseTwoIndex(self.nbasis, self.nfn)
        result._array[:] = (self._array+self._array.T)*factor/2.0
        return result

    def clear(self):
        '''Resets array to zeros element-wise.'''
        self._array[:] = 0.0

    def _iadd_array(self, arr, factor=1):
        '''Inplace addition of a numpy array. For matrix class internal use only
            ** Arguments: **
            arr
                An instance of a dense numpy ndarray
            factor
                A float or int for scaling
        '''
        if not isinstance(arr, np.ndarray):
            raise TypeError('The arr object must be a numpy nd.array instance. Got ', type(arr), ' instead.')
        self._array += arr*factor

    def iadd(self, other, factor=1):
        '''Inplace addition of other DenseTwoIndex
            ** Arguments: **
            other
                An instance of a DenseTwoIndex object
            factor
                A float or int for scaling
        '''
        self._array += other._array*factor

    def iaddnumber(self, number, factor=1):
        '''Inplace addition of number'''
        self._array[:] += number*factor

    def iaddt(self, other, factor=1):
        '''Inplace addition of transposed other DenseTwoIndex'''
        self._array += other._array.T*factor

    def iaddtosubmatrix(self, other, factor=1, *args):
        '''Inplace addition of other DenseTwoIndex'''
        ind = []
        for arg in args:
            ind.append(arg)
        self._array[ind[0]:ind[1], ind[2]:ind[3]] += other._array*factor

    def iaddouter(self, other, other2, factor=1):
        '''Inplace addition of outer product of two other DenseTwoIndex'''
        self._array += np.outer(other._array.ravel(), other2._array.ravel())*factor

    def iaddkron(self, other, other2, factor=1):
        '''Inplace addition of kronecker product of two other DenseTwoIndex'''
        self._array += np.kron(other._array, other2._array)*factor

    def contract(self, select, factor=1, *args):
        '''Contract any view of DenseTwoIndex'''
        ind = []
        for arg in args:
            ind.append(arg)
        if not ind:
            ind = [0, self._array.shape[0], 0, self._array.shape[1]]
        if select == 'aa':
            return np.einsum('aa', self._array[ind[0]:ind[1],ind[0]:ind[1]])*factor
        elif select == 'ab':
            return np.einsum('ab->...', self._array[ind[0]:ind[1],ind[2]:ind[3]])*factor
        else:
            raise NotImplementedError

    def contract_oneindex(self, other, other2, select, factor=1.0):
        '''Contract with other DenseTwoIndex and DenseOneIndex'''
        if select == 'b-ab':
            self._array[:] += np.einsum('ab,b->ab', other._array, other2._array.ravel())*factor
        elif select == 'a-ab':
            self._array[:] += np.einsum('ab,a->ab', other._array, other2._array.ravel())*factor
        else:
            raise NotImplementedError

    def contract_twoindex(self, other, factor=1.0):
        '''Contract with other DenseTwoIndex
            ** Arguments: **
                other
                    An instance of the DenseTwoIndex class
                factor
                    A float to scale the contraction by

        '''
        return (np.dot(self._array.ravel(), other._array.ravel()))*factor

    def contract_three2one(self, other3, other2, index, factor=1.0):
        ''' Add contraction of 3-index with 2-index tensor

            ** Arguments: **
                other3
                    An instance of the DenseThreeIndex class
                other2
                    An instance of the DenseTwoIndex class
                index
                    The index to contract over
                factor
                    A float to scale the contraction by

        '''
        if not isinstance(other3, DenseThreeIndex):
            raise TypeError('The other object must also be DenseThreeIndex instance. Got ', type(other3), ' instead.')

        if index == '13':
            self._array += other3.contract('ab->ac', other2._array)*factor
        elif index == '31':
            self._array += other3.contract('ab->ca', other2._array)*factor
        elif index == '21':
            self._array += other3.contract('bc->ba', other2._array)*factor
        elif index == '12':
            self._array += other3.contract('bc->ab', other2._array)*factor
        elif index == '132':
            self._array += other3.contract('cb->ac', other2._array)*factor
        else:
            raise NotImplementedError

    def iaddshift(self, lshift):
        '''Add positive shift to elements. If negative replace by shift'''
        pos = np.where(self._array >= 0)
        neg = np.where(self._array < 0)
        self._array[pos] += lshift
        self._array[neg] = lshift

    def expectation_value(self, dm, select=None):
        return np.dot(self._array.ravel(), dm._array.ravel())

    def trace(self):
        return np.trace(self._array)

    def trace_product(self, other):
        '''Return the trace of the product of this operator and of the given operator'''
        return np.dot(self._array.ravel(), other._array.T.ravel())

    def itranspose(self):
        '''In-place transpose'''
        self._array = self._array.T

    def iscale(self, factor):
        self._array *= factor

    def pass_array(self, ind=None):
        if ind:
            return self._array[ind]
        else:
            return self._array[:]

    def dot(self, vec0, vec1):
        return np.dot(vec0, np.dot(self._array, vec1))

    def dott(self, other):
        return np.dot(self._array, other._array.T)

    def idot(self, other):
        '''Inplace dot operator'''
        self._array[:] = np.dot(self._array, other._array)

    def iadddot(self, other, other2, factor=1.0):
        '''Inplace addition of dot product of two other DenseTwoIndex'''
        self._array[:] += np.dot(other._array, other2._array)*factor

    def iaddtdot(self, other, other2, factor=1.0):
        '''Inplace addition of dot product of one tansposed DenseTwoIndex'''
        self._array[:] += np.dot(other._array.T, other2._array)*factor

    def iadddott(self, other, other2, factor=1.0):
        '''Inplace addition of dot product of one tansposed DenseTwoIndex'''
        self._array[:] += np.dot(other._array, other2._array.T)*factor

    def multiply(self, other, vec, factor=1.0):
        '''Inplace addition of element-wise multiplication'''
        self._array[:] += other._array*vec._array*factor

    def multiplyt(self, other, vec, factor=1.0):
        '''Inplace addition of element-wise multiplication'''
        self._array[:] += other._array*vec._array.T*factor

    def multiplytt(self, other, vec, factor=1.0):
        '''Inplace addition of element-wise multiplication'''
        self._array[:] += other._array.T*vec._array.T*factor

    def imultiply(self, other, factor=1.0):
        '''Inplace element-wise multiplication'''
        self._array *= other._array*factor

    def imultiplyt(self, other, factor=1.0):
        '''Inplace transposed element-wise multiplication'''
        self._array *= other._array.T*factor

    def distance(self, other):
        '''Maximum difference between self and other two index object'''
        return abs(self._array.ravel() - other._array.ravel()).max()

    def iabsolute(self):
        self._array[:] = abs(self._array[:])

    def apply_basis_permutation(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        self._array[:] = self._array[permutation]
        self._array[:] = self._array[:,permutation]

    def apply_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.'''
        self._array *= signs
        self._array *= signs.reshape(-1,1)

    def apply_sorting_aabb(self, other):
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must be a DenseFourIndex instance. Got ', type(other), ' instead.')
        self._array[:] = other.get_slice('aabb->ab')
        return self._array

    def apply_sorting_abab(self, other):
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must be a DenseFourIndex instance. Got ', type(other), ' instead.')
        self._array[:] = other.get_slice('abab->ab')
        return self._array

    def apply_sorting_abba(self, other):
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must be a DenseFourIndex instance. Got ', type(other), ' instead.')
        self._array[:] = other.get_slice('abba->ab')
        return self._array

    def apply_sorting_pq(self, other):
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must be ac DenseFourIndex instance. Got ', type(other), ' instead.')
        self._array[:] = 2.0*other.get_slice('abab->ab') \
                                - other.get_slice('abba->ab')
        return self._array

    def apply_2index_trans(self, ao_integrals, aorb, aorb2=None):
        '''Perform two index transformation.
        '''
        if not aorb2:
            aorb2 = aorb
        self._array[:] = reduce(np.dot, [aorb.coeffs.T, ao_integrals._array, aorb2.coeffs])

    def compute_response_two_dm_ap1rog(self, one_dm, geminal, lagrange, select):
        '''Compute response 2-RDM for AP1roG

           ** Arguments **
           one_dm
               A 1DM. A TwoIndex instance.

           geminal
               The geminal coefficient matrix. A TwoIndex instance.

           lagrange
               The lagrange multipliers. A TwoIndex instance.

           select
               Either 'ppqq' or 'pqpq'. Note that the elements (iiii), i.e.,
               the 1DM, are stored in pqpq, while the elements (aaaa) are
               stored in ppqq.
        '''
        nocc = geminal._array.shape[0]
        ldotc = geminal.contract_twoindex(lagrange)
        if select == 'ppqq':
            self._array[:nocc,:nocc] = np.einsum("ja,ia->ij", lagrange._array, geminal._array)
            np.fill_diagonal(self._array[:nocc,:nocc],0)
            self._array[nocc:,nocc:] = np.einsum("ip,iq->pq", lagrange._array, geminal._array)
            self._array[nocc:,:nocc] = lagrange._array.T
            self._array[:nocc,nocc:] = geminal._array*(1-ldotc)
            self._array[:nocc,nocc:] += geminal._array*ldotc
            self._array[:nocc,nocc:] += 1*np.dot(geminal._array,np.dot(lagrange._array.T,geminal._array))
            self._array[:nocc,nocc:] -= 2*np.einsum("ip,ip->i",lagrange._array,geminal._array)[np.newaxis].T*geminal._array
            self._array[:nocc,nocc:] -= 2*np.einsum("ip,ip->p",lagrange._array,geminal._array)*geminal._array
            self._array[:nocc,nocc:] += 2*lagrange._array*geminal._array*geminal._array
        elif select == 'pqpq':
            for i in range(nocc):
                for j in range(i+1,nocc):
                    tmp = np.dot(lagrange._array[i,:],geminal._array[i,:].T)+np.dot(lagrange._array[j,:],geminal._array[j,:].T)
                    self._array[i,j] = 1.0-tmp
                    self._array[j,i] = 1.0-tmp
                self._array[i,i] = 1.0-np.dot(lagrange._array[i,:],geminal._array[i,:].T)
            self._array[:nocc,nocc:] = one_dm._array[nocc:].T-lagrange._array*geminal._array
            self._array[nocc:,:nocc] = one_dm._array[nocc:]-geminal._array.T*lagrange._array.T
        else:
            raise NotImplementedError

    def compute_ps2_two_dm_ap1rog(self, one_dm, geminal, lagrange, select):
        '''Compute PS2-type 2-RDM for AP1roG.
        '''
        nocc = geminal._array.shape[0]
        ldotc = geminal.contract_twoindex(lagrange)
        value = 1+ldotc
        if select == 'ppqq':
            self._array[:nocc,:nocc] = np.einsum("ia,ja->ij", lagrange._array, geminal._array)
            np.fill_diagonal(self._array[:nocc,:nocc],0)
            self._array[nocc:,nocc:] = np.einsum("ip,iq->pq", lagrange._array, geminal._array)
            self._array[nocc:,:nocc] = lagrange._array.T
            self._array[:nocc,nocc:] = geminal._array
            self._array[:nocc,nocc:] += geminal._array*ldotc
            self._array[:nocc,nocc:] += 1*np.dot(geminal._array,np.dot(lagrange._array.T,geminal._array))
            self._array[:nocc,nocc:] -= 2*np.einsum("ip,ip->i",lagrange._array,geminal._array)[np.newaxis].T*geminal._array
            self._array[:nocc,nocc:] -= 2*np.einsum("ip,ip->p",lagrange._array,geminal._array)*geminal._array
            self._array[:nocc,nocc:] += 2*lagrange._array*geminal._array*geminal._array
        elif select == 'pqpq':
            for i in range(nocc):
                for j in range(i+1,nocc):
                    tmp = np.dot(lagrange._array[i,:],geminal._array[i,:].T)+np.dot(lagrange._array[j,:],geminal._array[j,:].T)
                    self._array[i,j] = value-tmp
                    self._array[j,i] = value-tmp
                self._array[i,i] = value-np.dot(lagrange._array[i,:],geminal._array[i,:].T)
            self._array[:nocc,nocc:] = one_dm._array[nocc:].T-lagrange._array*geminal._array
            self._array[nocc:,:nocc] = one_dm._array[nocc:]-geminal._array.T*lagrange._array.T
        else:
            raise NotImplementedError

    def compute_response_four_dm_ap1rog(self, two_dm, geminal, lagrange, select):
        '''Compute response 4-RDM for AP1roG'''
        nocc = geminal._array.shape[0]
        if select == 'udud':
            self._array[:nocc,:nocc] = two_dm._array[:nocc, :nocc]
            self._array[:nocc,nocc:] = two_dm._array[:nocc, nocc:]
            self._array[nocc:,:nocc] = two_dm._array[:nocc, nocc:].T
            np.fill_diagonal(self._array[:nocc,:nocc], 0)


class DenseThreeIndex(ThreeIndex):
    """Dense three-dimensional object.

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
        self._array = np.zeros((nbasis, nbasis, nbasis), float)
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            log.mem.denounce(self._array.nbytes)

    def __check_init_args__(self, nbasis):
        assert nbasis == self.nbasis

    @classmethod
    def from_hdf5(cls, grp):
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

    def set_element(self, i, j, k, value):
        self._array[i,j,k] = value

    def get_element(self, i, j, k):
        return self._array[i,j, k]

    def clear(self):
        self._array[:] = 0.0

    def iscale(self, factor):
        '''Inplace scaling by factor'''
        self._array[:] = self._array[:]*factor

    def assign_array(self, other):
        ''''''
        self._array[:] = other[:]

    def contract(self, indices, arr):
        ''' Contracts a dense 2-dim array with this DenseThreeIndex term.

            ** Arguments: **

            indices
                A string with two sets of letters and "->" separating them.
                These are the indices of the 2-dim array and the output.
                The indices of this DenseThreeIndex term will be "abc".
                See numpy einsum notation documentation for more details.

            arr
                The 2-dim dense numpy array to be contracted over.

        '''
        indc_in,indc_out = "".join(indices.split()).split("->")
        assert len(indc_in) == 2 and len(indc_out) == 2
        #TODO: add safety checks for indices? Einsum already checks...

        return np.einsum('abc,'+indices, self._array,arr)

    def expand_tothreeindex(self, other, other2, select, factor=1.0):
        '''Expand One and OneIndex to ThreeIndex.

           **Arguments:**

           other
                A DenseTwoIndex object.

           other2
                A DenseOneIndex object.

           select:
                Contraction type
        '''
        if select == 'cab':
            self._array[:] += np.einsum('ab,c->cab', other._array, other2._array.ravel())*factor
        elif select == 'acb':
            self._array[:] += np.einsum('ab,c->acb', other._array, other2._array.ravel())*factor
        else:
            raise NotImplementedError

    def contract_twoindex(self, other, other2, select, factor=1.0):
        '''Contract other ThreeIndex with TwoIndex.

           **Arguments:**

           other
                A DenseThreeIndex object.

           other2
                A DenseTwoIndex object.

           select:
                Contraction type
        '''
        if select == 'db-adc':
            self._array[:] += np.einsum('abc,db->adc', other._array, other2._array)*factor
        elif select == 'dc-adb':
            self._array[:] += np.einsum('abc,dc->adb', other._array, other2._array)*factor
        elif select == 'db-dac':
            self._array[:] += np.einsum('abc,db->dac', other._array, other2._array)*factor
        elif select == 'dc-dab':
            self._array[:] += np.einsum('abc,dc->dab', other._array, other2._array)*factor
        else:
            raise NotImplementedError

    def apply_sorting_ijkk(self, other, switch=False):
        ''' Sort four-index integrals for OAP1roG (<ij|kk> -> <ijk>)

            ** Arguments: **

            other
                An instance of a DenseFourIndex object.
            switch
                A boolean #FIXME

        '''
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must also be DenseFourIndex instance. Got ', type(other), ' instead.')
        if switch:
            other.get_slice('abcc->bac', out=self._array)
        else:
            other.get_slice('abcc->abc', out=self._array)

    def apply_sorting_pqrq_woexchange(self, other):
        '''Sort four-index integrals for OAP1roG (<ij|kj> -> <ijk>)
        '''
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must also be DenseFourIndex instance.')
        other.get_slice('abcb->abc', out=self._array)

    def apply_sorting_pqrq_exchange(self, other):
        '''Sort four-index exchange integrals for OAP1roG (<ij||kj> -> <ijk>)
        '''
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must also be DenseFourIndex instance.')
        other.get_slice('abcb->abc', out=self._array)
        other.get_slice('abbc->abc', out=self._array, subtract=True, clear=False)
        return self._array

    def apply_sorting_pqrq(self, other):
        '''Sort four-index exchange integrals for OAP1roG (<ij||kj> -> <ijk>)
        '''
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must also be DenseFourIndex instance. Got ', type(other), ' instead.')
        other.get_slice('abcb->abc', out=self._array)
        self._array *= 2.0
        other.get_slice('abbc->abc', out=self._array, subtract=True, clear=False)
        return self._array

    def compute_response_three_dm_ap1rog(self, geminal, lagrange, select):
        '''Compute response 3-RDM for AP1roG.
           Only the symmetry-unique elements are calculated.
           Only one block is stored.

           **Arguments:**

           geminal
                A DenseTwoIndex Geminal coefficient matrix.

           lagrange
                A DenseTwoIndex instance of Lagrange multipliers.

           select
                Block of 3-DM; uuu, uud, uudoff, udd, uddoff
        '''
        nocc = geminal._array.shape[0]
        nvirt = geminal._array.shape[1]
        ldotc = geminal.contract_twoindex(lagrange)
        nd_lagrange = lagrange._get_dense()
        nd_geminal = geminal._get_dense()
        if select == 'uuu':
            unity = np.ones(nocc)
            value = 1.0
            tmpo = np.dot(nd_lagrange, nd_geminal.T).diagonal()
            tmpv = np.dot(nd_lagrange.T, nd_geminal).diagonal()
            # ooo-block
            self._array[:nocc, :nocc, :nocc] = value
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", tmpo, unity, unity)
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", unity, tmpo, unity)
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", unity, unity, tmpo)
            # voo-block
            self._array[nocc:, :nocc, :nocc] = np.einsum("a...,b...,c...->abc", tmpv, unity, unity)
            self._array[nocc:, :nocc, :nocc] -= np.einsum("ab...,c...->abc", (nd_lagrange.T*nd_geminal.T), unity)
            self._array[nocc:, :nocc, :nocc] -= np.einsum("ac...,b...->abc", (nd_lagrange.T*nd_geminal.T), unity)
            # ovo-block
            self._array[:nocc, nocc:, :nocc] = np.einsum("abc->bac", self._array[nocc:, :nocc, :nocc])
            # oov-block
            self._array[:nocc, :nocc, nocc:] = np.einsum("abc->bca", self._array[nocc:, :nocc, :nocc])
            # Set [i,i,k] equal to zero and account for double-subtraction
            # This is not a nice solution:
            # TODO: find better way
            for i in range(nocc):
                for j in range(nocc):
                    self._array[i,i,j] = 0.0
                    self._array[i,j,i] = 0.0
                    self._array[j,i,i] = 0.0
                self._array[i,i,i] = 0.0
                for a in range(self.nbasis-nocc):
                    self._array[i,i,(a+nocc)] = 0.0
                    self._array[i,(a+nocc),i] = 0.0
                    self._array[(a+nocc),i,i] = 0.0
        elif select == 'uud':
            unity = np.ones(nocc)
            value = 1.0
            tmpo = np.dot(nd_lagrange, nd_geminal.T).diagonal()
            tmpv = np.dot(nd_lagrange.T, nd_geminal).diagonal()
            # ooo-block
            self._array[:nocc, :nocc, :nocc] = value
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", tmpo, unity, unity)
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", unity, tmpo, unity)
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", unity, unity, tmpo)
            # voo-block
            self._array[nocc:, :nocc, :nocc] = np.einsum("a...,b...,c...->abc", tmpv, unity, unity)
            self._array[nocc:, :nocc, :nocc] -= np.einsum("ab...,c...->abc", (nd_lagrange.T*nd_geminal.T), unity)
            self._array[nocc:, :nocc, :nocc] -= np.einsum("ac...,b...->abc", (nd_lagrange.T*nd_geminal.T), unity)
            # ovo-block
            self._array[:nocc, nocc:, :nocc] = np.einsum("abc->bac", self._array[nocc:, :nocc, :nocc])
            # oov-block
            self._array[:nocc, :nocc, nocc:] = np.einsum("abc->cba", self._array[nocc:, :nocc, :nocc])
            # Set [i,i,k] equal to zero and account for double-subtraction
            # This is not a nice solution:
            # TODO: find better way
            for i in range(nocc):
                for j in range(nocc):
                    self._array[i,j,j] += tmpo[j]
                    self._array[j,i,j] += tmpo[j]
                    self._array[i,i,j] = 0.0
                for a in range(self.nbasis-nocc):
                    self._array[(a+nocc),i,i] += nd_lagrange[i,a]*nd_geminal[i,a]
                    self._array[i,(a+nocc),i] += nd_lagrange[i,a]*nd_geminal[i,a]
                    self._array[i,i,(a+nocc)] = 0.0
                    # ovv-block
                    self._array[i,(a+nocc),(a+nocc)] = tmpv[a]-nd_lagrange[i,a]*nd_geminal[i,a]
                    # vov-block
                    self._array[(a+nocc),i,(a+nocc)] = tmpv[a]-nd_lagrange[i,a]*nd_geminal[i,a]
        elif select == 'uudoff':
            unity = np.ones(nocc)
            tmp = np.zeros((nocc,nvirt))
            tmpv = np.dot(nd_geminal, nd_lagrange.T)
            tmpo = np.dot(nd_lagrange.T, nd_geminal)
            # ooo-block
            self._array[:nocc, :nocc, :nocc] = np.einsum("a...,bc...->abc", unity, tmpv)
            # voo-block
            self._array[nocc:, :nocc, :nocc] = np.einsum("ab...,cb...->bac", nd_geminal, nd_lagrange)
            # ovo-block
            self._array[:nocc, nocc:, :nocc] = np.einsum("a...,bc...->acb", unity, nd_lagrange)
            # oov-block
            tmp += nd_geminal*(1-ldotc)
            tmp += nd_geminal*ldotc
            tmp += 1*np.dot(nd_geminal,np.dot(nd_lagrange.T,nd_geminal))
            tmp -= 2*np.einsum("ip,ip->i",nd_lagrange,nd_geminal)[np.newaxis].T*nd_geminal
            tmp -= 2*np.einsum("ip,ip->p",nd_lagrange,nd_geminal)*nd_geminal
            tmp += 2*nd_lagrange*nd_geminal*nd_geminal
            self._array[:nocc, :nocc, nocc:] = np.einsum("a...,bc...->abc", unity, tmp)
            # ovv-block
            self._array[:nocc, nocc:, nocc:] = np.einsum("a...,bc...->abc", unity, tmpo)
            self._array[:nocc, nocc:, nocc:] -= np.einsum("ab...,ac...->abc", nd_lagrange, nd_geminal)
            # vov-block
            # Set [i,i,k] equal to zero and account for double-subtraction
            # This is not a nice solution:
            # TODO: find better way
            for i in range(nocc):
                for j in range(nocc):
                    self._array[i,j,j] = 0.0
                    self._array[j,i,j] = 0.0
                    self._array[i,i,j] = 0.0
                self._array[i,i,i] = 0.0
                for a in range(self.nbasis-nocc):
                    self._array[i,i,(a+nocc)] = 0.0
                    # ovo-block
                    self._array[i,(a+nocc),i] = 0.0
                    # voo-block
                    self._array[(a+nocc),i,i] = 0.0
        elif select == 'udd':
            unity = np.ones(nocc)
            value = 1.0
            tmpo = np.dot(nd_lagrange, nd_geminal.T).diagonal()
            tmpv = np.dot(nd_lagrange.T, nd_geminal).diagonal()
            # ooo-block
            self._array[:nocc, :nocc, :nocc] = value
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", tmpo, unity, unity)
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", unity, tmpo, unity)
            self._array[:nocc, :nocc, :nocc] -= np.einsum("a...,b...,c...->abc", unity, unity, tmpo)
            # voo-block
            self._array[nocc:, :nocc, :nocc] = np.einsum("a...,b...,c...->abc", tmpv, unity, unity)
            self._array[nocc:, :nocc, :nocc] -= np.einsum("ab...,c...->abc", (nd_lagrange.T*nd_geminal.T), unity)
            self._array[nocc:, :nocc, :nocc] -= np.einsum("ac...,b...->abc", (nd_lagrange.T*nd_geminal.T), unity)
            # ovo-block
            self._array[:nocc, nocc:, :nocc] = np.einsum("a...,b...,c...->abc", unity, tmpv, unity)
            self._array[:nocc, nocc:, :nocc] -= np.einsum("ab...,c...->abc", (nd_lagrange*nd_geminal), unity)
            self._array[:nocc, nocc:, :nocc] -= np.einsum("a...,bc...->abc", unity, (nd_lagrange.T*nd_geminal.T))
            # oov-block
            self._array[:nocc, :nocc, nocc:] = np.einsum("abc->acb", self._array[:nocc, nocc:, :nocc])
            # Set [i,i,k] equal to zero and account for double-subtraction
            # This is not a nice solution:
            # TODO: find better way
            for i in range(nocc):
                for j in range(nocc):
                    self._array[j,j,i] += tmpo[j]
                    self._array[j,i,j] += tmpo[j]
                    self._array[i,j,j] = 0.0
                for a in range(self.nbasis-nocc):
                    self._array[i,(a+nocc),i] += nd_lagrange[i,a]*nd_geminal[i,a]
                    self._array[i,i,(a+nocc)] += nd_lagrange[i,a]*nd_geminal[i,a]
                    self._array[(a+nocc),i,i] = 0.0
                    # ovv-block
                    self._array[i,(a+nocc),(a+nocc)] = tmpv[a]-nd_lagrange[i,a]*nd_geminal[i,a]
                    # vov-block
                    self._array[(a+nocc),i,(a+nocc)] = tmpv[a]-nd_lagrange[i,a]*nd_geminal[i,a]
        elif select == 'uddoff':
            unity = np.ones(nocc)
            tmp = np.zeros((nocc,nvirt))
            tmpv = np.dot(nd_geminal, nd_lagrange.T)
            tmpo = np.dot(nd_lagrange.T, nd_geminal)
            # ooo-block
            self._array[:nocc, :nocc, :nocc] = np.einsum("ab...,c...->abc", tmpv, unity)
            # voo-block
            self._array[nocc:, :nocc, :nocc] = np.einsum("ab...,c...->bac", nd_lagrange, unity)
            # ovo-block
            tmp += nd_geminal*(1-ldotc)
            tmp += nd_geminal*ldotc
            tmp += 1*np.dot(nd_geminal,np.dot(nd_lagrange.T,nd_geminal))
            tmp -= 2*np.einsum("ip,ip->i",nd_lagrange,nd_geminal)[np.newaxis].T*nd_geminal
            tmp -= 2*np.einsum("ip,ip->p",nd_lagrange,nd_geminal)*nd_geminal
            tmp += 2*nd_lagrange*nd_geminal*nd_geminal
            self._array[:nocc, nocc:, :nocc] = np.einsum("ab...,c...->abc", tmp, unity)
            # vvo-block
            self._array[nocc:, nocc:, :nocc] = np.einsum("ab...,c...->abc", tmpo, unity)
            self._array[nocc:, nocc:, :nocc] -= np.einsum("ac...,bc...->abc", nd_lagrange.T, nd_geminal.T)
#           # Set [i,i,k] equal to zero and account for double-subtraction
#           # This is not a nice solution:
#           # TODO: find better way
            for i in range(nocc):
                for j in range(nocc):
                    self._array[i,j,j] = 0.0
                    self._array[j,i,j] = 0.0
                    self._array[i,i,j] = 0.0
                    self._array[i,i,i] = 0.0
                for a in range(self.nbasis-nocc):
                    # ovo-block
                    self._array[i,(a+nocc),i] = 0.0
                    # voo-block
                    self._array[(a+nocc),i,i] = 0.0
        else:
            raise NotImplementedError


class DenseFourIndex(FourIndex):
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
        self._array = np.zeros((nbasis, nbasis, nbasis, nbasis))
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)

    def __check_init_args__(self, nbasis):
        assert nbasis == self.nbasis

    @classmethod
    def from_hdf5(cls, grp):
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
        return self._array[i, j, k, l]

    def get_slice(self, indices, out=None, subtract=False, clear=False):
        """Returns a numpy array 2-index slice of the fourindex object.

            ** Arguments **
                indc_in
                    A string of length 4, comprised of some combination of letters.
                    These are the indices that will be read from the two electron integrals.
                    The letters "x,y,z" are reserved for internal use.
                    See numpy einstein summation documentation for an example.
                indc_out
                    A string, comprised of the same letters in indc_in.
                    These are the output indices. Reversed order will give a transpose.
                    See numpy einstein summation documentation for an example.
        """
        #error checking
        indc_in,indc_out = "".join(indices.split()).split("->")
        assert len(indc_in) == 4
        assert set(indc_in) == set(indc_out)
        assert len(set("xyz") & set(indc_in)) == 0
        if out is not None:
            if clear:
                out[:] = 0.0
            if subtract:
                out -= np.einsum(indices,self._array)
            else:
                out += np.einsum(indices, self._array)
            return out
        else:
            return np.einsum(indices, self._array)

    def assign(self, other):
        if not isinstance(other, DenseFourIndex):
            raise TypeError('The other object must also be DenseFourIndex instance. Got ', type(other), ' instead.')
        self._array[:] = other._array

    def assign_array(self, other):
        '''Assign a dense numpy array to the underlying representation.'''
        self._array[:] = other[:]

    def copy(self):
        '''Return a copy of the current four-index operator'''
        result = DenseFourIndex(self.nbasis)
        result.assign(self)
        return result

    def iscale(self, factor):
        self._array *= factor

    def check_symmetry(self):
        """Check the symmetry of the array."""
        assert abs(self._array - self._array.transpose(1,0,3,2)).max() == 0.0
        assert abs(self._array - self._array.transpose(2,3,0,1)).max() == 0.0
        assert abs(self._array - self._array.transpose(3,2,1,0)).max() == 0.0
        assert abs(self._array - self._array.transpose(2,1,0,3)).max() == 0.0
        assert abs(self._array - self._array.transpose(3,0,1,2)).max() == 0.0
        assert abs(self._array - self._array.transpose(0,3,2,1)).max() == 0.0
        assert abs(self._array - self._array.transpose(1,2,3,0)).max() == 0.0

    def transpose(self):
        self._array = self._array.transpose(1,0,3,2)

    def esum(self):
        return np.sum(self._array)

    def assign_hessian(self, other1, other2, other3, other4, other5):
        for i in range(self.nbasis):
            self._array[:,i,i,:] += other1._array
            self._array[i,:,:,i] += other1._array
            self._array[i,:,i,:] -= other1._array
            self._array[:,i,:,i] -= other1._array
            self._array[i,:,i,:] -= other2._array[i,:,:]
            self._array[i,:,:,i] += other2._array[i,:,:]
            self._array[:,i,:,i] -= other3._array[:,i,:]
            self._array[:,i,i,:] += other3._array[:,i,:]
            self._array[:,i,i,:] += other4._array[:,i,:]
            self._array[:,i,:,i] -= other4._array[:,i,:]
            self._array[i,:,:,i] += other5._array[i,:,:]
            self._array[i,:,i,:] -= other5._array[i,:,:]

    def contract_twoindex(self, other, other2, select, factor=1.0):
        '''Contract other FourIndex with TwoIndex.

           **Arguments:**

           other
                A DenseFourIndex object.

           other2
                A DenseTwoIndex object.

           select:
                Contraction type
        '''
        if select == 'cd-acbd':
            self._array[:] += other.contract('cd->acbd', other2._array)*factor
        elif select == 'cd-acdb':
            self._array[:] += other.contract('cd->acdb', other2._array)*factor
        elif select == 'cb-acdb':
            self._array[:] += other.contract('cb->acdb', other2._array)*factor
        elif select == 'cb-acbd':
            self._array[:] += other.contract('cb->acbd', other2._array)*factor
        elif select == 'ab-acbd':
            self._array[:] += other.contract('ab->acbd', other2._array)*factor
        elif select == 'ab-acdb':
            self._array[:] += other.contract('ab->acdb', other2._array)*factor
        elif select == 'ad-acbd':
            self._array[:] += other.contract('ad->acbd', other2._array)*factor
        elif select == 'ad-acdb':
            self._array[:] += other.contract('ad->acdb', other2._array)*factor
        elif select == 'ad-abcd':
            self._array[:] += other.contract('ad->abcd', other2._array)*factor
        elif select == 'ad-abdc':
            self._array[:] += other.contract('ad->abdc', other2._array)*factor
        elif select == 'bd-abcd':
            self._array[:] += other.contract('bd->abcd', other2._array)*factor
        elif select == 'bd-abdc':
            self._array[:] += other.contract('bd->abdc', other2._array)*factor
        elif select == 'bc-abdc':
            self._array[:] += other.contract('bc->abdc', other2._array)*factor
        elif select == 'bc-abcd':
            self._array[:] += other.contract('bc->abcd', other2._array)*factor
        elif select == 'ac-abcd':
            self._array[:] += other.contract('ac->abcd', other2._array)*factor
        elif select == 'ac-abdc':
            self._array[:] += other.contract('ac->abdc', other2._array)*factor
        else:
            raise NotImplementedError

    def contract(self, indices, arr):
        ''' Contracts a dense 2-dim array with this DenseTwoIndex term.

            ** Arguments: **

            indices
                A string with two sets of letters and "->" separating them.
                These are the indices of the 2-dim array and the output.
                The indices of this DenseFourIndex term will be "abcd".
                See numpy einsum notation documentation for more details.

            arr
                The 2-dim dense numpy array to be contracted over.

        '''
        indc_in,indc_out = "".join(indices.split()).split("->")
        assert len(indc_in) == 2 and len(indc_out) == 4
        #TODO: add safety checks for indices? Einsum already checks...

        return np.einsum('abcd,'+indices, self._array, arr)

    def apply_direct(self, dm, output):
        """Compute the direct dot product with a density matrix."""
        if not isinstance(dm, DenseTwoIndex):
            raise TypeError('The dm argument must be a DenseTwoIndex class. Got ', type(dm), ' instead.')
        if not isinstance(output, DenseTwoIndex):
            raise TypeError('The output argument must be a DenseTwoIndex class. Got ', type(output), ' instead.')
        output.assign_array(np.tensordot(self._array, dm._get_dense(), ([1,3], [1,0])))

    def apply_exchange(self, dm, output):
        """Compute the exchange dot product with a density matrix."""
        if not isinstance(dm, DenseTwoIndex):
            raise TypeError('The dm argument must be a DenseTwoIndex class. Got ', type(dm), ' instead.')
        if not isinstance(output, DenseTwoIndex):
            raise TypeError('The output argument must be a DenseTwoIndex class. Got ', type(output), ' instead.')
        output.assign_array(np.tensordot(self._array, dm._get_dense(), ([1,2], [1,0])))

    def add_exchange_part(self):
        '''Sort four-index exchange integrals for OAP1roG (<ij||kj> -> <ijk>)
        '''
        tmp = np.einsum('abcd->abdc', self._array)
        self._array[:] = self._array-tmp
        del tmp
        return self._array

    def apply_four_index_transform_tensordot(self, ao_integrals, aorb, aorb2=None, aorb3=None, aorb4=None):
        '''Perform four index transformation using np.tensordot.
        '''
        if not aorb2:
            aorb2 = aorb
        if not aorb3:
            aorb3 = aorb
        if not aorb4:
            aorb4 = aorb

        if not isinstance(ao_integrals, DenseFourIndex):
            raise TypeError('The AO integral argument must be a DenseFourIndex class. Got ', type(ao_integrals), ' instead.')
        self._array[:] = np.tensordot(ao_integrals._array,aorb2.coeffs,axes=([0],[0]))
        self._array[:] = np.tensordot(self._array,aorb.coeffs,axes=([0],[0]))
        self._array[:] = np.tensordot(self._array,aorb4.coeffs,axes=([0],[0]))
        self._array[:] = np.tensordot(self._array,aorb3.coeffs,axes=([0],[0]))

    def apply_four_index_transform_einsum(self, ao_integrals, aorb, aorb2=None, aorb3=None, aorb4=None):
        '''Perform four index transformation using np.einsum.
        '''
        if not aorb2:
            aorb2 = aorb
        if not aorb3:
            aorb3 = aorb
        if not aorb4:
            aorb4 = aorb

        if not isinstance(ao_integrals, DenseFourIndex):
            raise TypeError('The AO integral argument must be a DenseFourIndex class. Got ', type(ao_integrals), ' instead.')
        self._array[:] = np.einsum('st,pqrt->pqrs',aorb.coeffs.T,ao_integrals._array,order='C',casting='no')
        self._array[:] = np.einsum('rs,pqst->pqrt',aorb2.coeffs.T,self._array,casting='no',order='C')
        self._array[:] = np.einsum('qr,prst->pqst',aorb3.coeffs.T,self._array,casting='no',order='C')
        self._array[:] = np.einsum('ab,bqrs->aqrs',aorb4.coeffs.T,self._array,casting='no',order='C')

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
