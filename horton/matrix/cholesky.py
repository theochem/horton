# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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
"""Cholesky decomposition of four-index objects"""


import numpy as np

from horton.utils import check_type, check_options, doc_inherit
from horton.matrix.base import parse_four_index_transform_exps, FourIndex
from horton.matrix.cext import slice_to_three_abbc_abc, \
    slice_to_three_abcc_bac, slice_to_three_abcc_abc
from horton.matrix.dense import DenseLinalgFactory, DenseExpansion, \
    DenseTwoIndex, DenseThreeIndex, DenseFourIndex


__all__ = [
    'CholeskyFourIndex', 'CholeskyLinalgFactory',
]


class CholeskyLinalgFactory(DenseLinalgFactory):
    @doc_inherit(DenseLinalgFactory)
    def create_four_index(self, nbasis=None, nvec=None, array=None, array2=None):
        nbasis = nbasis or self.default_nbasis
        return CholeskyFourIndex(nbasis, nvec, array, array2)

    @doc_inherit(DenseLinalgFactory)
    def _check_four_index_init_args(self, four_index, nbasis=None, nvec=None, array=None):
        nbasis = nbasis or self.default_nbasis
        four_index.__check_init_args__(nbasis, nvec)

    create_four_index.__check_init_args__ = _check_four_index_init_args


class CholeskyFourIndex(FourIndex):
    """Cholesky symmetric four-dimensional matrix.
    """

    #
    # Constructor and destructor
    #

    def __init__(self, nbasis, nvec=None, array=None, array2=None):
        """
           **Arguments:**

           nbasis
                The number of basis functions.

           **Optional arguments:**

           nvec
                The number of (2-index) Cholesky vectors.

           array
                The array with Cholesky vectors, shape = (nvec, nbasis, nbasis).

           array2
                The second set of Cholesky vectors, if different from the first.

           Either nvec or array must be given (or both).
        """
        def check_array(a, name):
            if a.ndim != 3:
                raise TypeError('Argument %s has %i dimensions, expecting 3.' % (name, a.ndim))
            if nvec is not None and nvec != a.shape[0]:
                raise TypeError('nvec does not match %s.shape[0].' % name)
            if not (nbasis == a.shape[1] and nbasis == a.shape[2]):
                raise TypeError('nbasis does not match %s.shape[1] or %s.shape[2].' % (name, name))

        if array is None:
            self._self_alloc = True
            if nvec is None:
                raise TypeError('Either nvec or array must be given (or both).')
            if array2 is not None:
                raise TypeError('Argument array2 only allowed when array is given.')
            self._array = np.zeros([nvec, nbasis, nbasis])
            self._array2 = self._array
        else:
            self._self_alloc = False
            check_array(array, 'array')
            self._array = array
            if array2 is None:
                self._array2 = self._array
            else:
                check_array(array2, 'array2')
                self._array2 = array2

    #
    # Properties
    #

    def _get_shape(self):
        '''The shape of the object'''
        return (self._array.shape[1], self._array2.shape[1], self._array.shape[2], self._array2.shape[2])

    shape = property(_get_shape)

    #
    # Methods from base class
    #

    def __check_init_args__(self, nbasis, nvec):
        '''Is self compatible with the given constructor arguments?'''
        assert self._array is not None
        assert nbasis == self.nbasis
        assert nvec == self.nvec

    def __eq__(self, other):
        '''Compare self with other'''
        return isinstance(other, CholeskyFourIndex) and \
            other.nbasis == self.nbasis and \
            other.nvec == self.nvec and \
            other.is_decoupled == self.is_decoupled and \
            (other._array == self._array).all() and \
            (other._array2 == self._array2).all()

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct an instance from data previously stored in an h5py.Group.

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        nvec = grp['array'].shape[0]
        nbasis = grp['array'].shape[1]
        result = cls(nbasis, nvec)
        grp['array'].read_direct(result._array)
        if 'array2' in grp:
            result.decouple_array2()
            grp['array2'].read_direct(result._array2)
        return result

    def to_hdf5(self, grp):
        '''Dump this object in an h5py.Group

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        grp.attrs['class'] = self.__class__.__name__
        grp['array'] = self._array
        if self._array is not self._array:
            grp['array2'] = self._array2

    def new(self):
        '''Return a new four-index object with the same nbasis'''
        return CholeskyFourIndex(self.nbasis, self.nvec)

    def _check_new_init_args(self, other):
        '''Check whether an already initialized object is compatible'''
        other.__check_init_args__(self.nbasis, self.nvec)

    new.__check_init_args__ = _check_new_init_args

    def clear(self):
        '''Reset all elements to zero.'''
        self._array[:] = 0.0
        if self._array is not self._array2:
            self._array2[:] = 0.0

    def copy(self):
        '''Return a copy of the current four-index operator'''
        result = CholeskyFourIndex(self.nbasis, self.nvec)
        result.assign(self)
        return result

    def assign(self, other):
        '''Assign with the contents of another object

           **Arguments:**

           other
                Another CholeskyFourIndex object.
        '''
        check_type('other', other, CholeskyFourIndex)
        self._array[:] = other._array
        if other._array is other._array2:
            self.reset_array2()
        else:
            self.decouple_array2()
            self._array2[:] = other._array2

    def randomize(self):
        '''Fill with random normal data'''
        self._array[:] = np.random.normal(0, 1, self._array.shape)
        if self.is_decoupled:
            self._array2[:] = np.random.normal(0, 1, self._array2.shape)

    def permute_basis(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        # Easy enough but irrelevant
        raise NotImplementedError

    def change_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.'''
        # Easy enough but irrelevant
        raise NotImplementedError

    def iadd(self, other, factor):
        '''This method is not supported due to the Cholesky decomposition.'''
        raise NotImplementedError

    def iscale(self, factor):
        '''In-place multiplication with a scalar

           **Arguments:**

           factor
                A scalar factor.
        '''
        self._array *= np.sqrt(factor)

        if self._array is not self._array2:
            #arrays have been transformed
            self._array2 *= np.sqrt(factor)

    def get_element(self, i, j, k, l):
        '''Return a matrix element'''
        return np.dot(self._array[:,i,k], self._array2[:,j,l])

    def set_element(self, i, j, k, l, value):
        '''This method is not supported due to the Cholesky decomposition.'''
        raise NotImplementedError

    #
    # Properties
    #

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._array.shape[1]

    nbasis = property(_get_nbasis)

    def _get_nvec(self):
        '''The number of Cholesky vectors'''
        return self._array.shape[0]

    nvec = property(_get_nvec)

    def _get_is_decoupled(self):
        return self._array is not self._array2

    is_decoupled = property(_get_is_decoupled)

    #
    # New methods for this implementation
    # TODO: consider adding these to base class
    #

    def decouple_array2(self):
        '''Allocates a second Cholesky vector if not done yet'''
        if self._array2 is self._array:
            self._array2 = self._array.copy()

    def reset_array2(self):
        """Deallocates the second cholesky vector and sets it to match the first.
        """
        if self._array2 is not self._array:
            self._array2 = self._array

    def get_dense(self):
        '''Return the DenseFourIndex equivalent. ONLY FOR TESTING. SUPER SLOW.
        '''
        result = DenseFourIndex(self.nbasis)
        np.einsum('kac,kbd->abcd', self._array, self._array2, out=result._array)
        return result

    def is_symmetric(self, symmetry=2, rtol=1e-5, atol=1e-8):
        '''Check the symmetry of the array.

           **Optional arguments:**

           symmetry
                The symmetry to check. See :ref:`dense_matrix_symmetry`
                for more details.

           rtol and atol
                relative and absolute tolerance. See to ``np.allclose``.
        '''
        if self.is_decoupled and symmetry in (2, 8):
            return False
        if symmetry in (4, 8):
            if not np.allclose(self._array, self._array.swapaxes(1,2), rtol, atol):
                return False
            if self.is_decoupled and not np.allclose(self._array2, self._array2.swapaxes(1,2), rtol, atol):
                return False
        return True

    def symmetrize(self, symmetry=8):
        check_options('symmetry', symmetry, 1, 2, 4, 8)
        if symmetry in (2, 8) and self.is_decoupled:
            # This is a different type of symmetrization than in the dense case!
            self._array[:] += self._array2
            self._array *= 0.5
            self.reset_array2()
        if symmetry in (4, 8):
            self._array[:] = self._array + self._array.transpose(0,2,1)
            if self.is_decoupled:
                self._array2[:] = self._array2 + self._array2.transpose(0,2,1)

    def itranspose(self):
        '''In-place transpose: ``0,1,2,3 -> 1,0,3,2``'''
        if self.is_decoupled:
            self._array, self._array2 = self._array2, self._array

    def sum(self):
        '''Return the sum of all elements. EXPENSIVE!'''
        return np.tensordot(self._array, self._array2,(0,0)).sum() #expensive!!

    def iadd_exchange(self):
        '''In-place addition of its own exchange contribution'''
        raise NotImplementedError

    def slice_to_two(self, subscripts, out=None, factor=1.0, clear=True):
        """Returns a two-index contraction of the four-index object.

           **Arguments:**

           subscripts
                Any of ``aabb->ab``, ``abab->ab``, ``abba->ab``

           **Optional arguments:**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`
        """
        # Error checking
        check_options('subscripts', subscripts, 'aabb->ab', 'abab->ab', 'abba->ab')
        # Handle output argument
        if out is None:
            out = DenseTwoIndex(self.nbasis)
        else:
            check_type('out', out, DenseTwoIndex)
            if clear:
                out.clear()
        # Actual computation
        if subscripts == 'aabb->ab':
            out._array[:] += factor*np.einsum('xab,xab->ab', self._array, self._array2)
        elif subscripts == 'abab->ab':
            out._array[:] += factor*np.einsum('xaa,xbb->ab', self._array, self._array2)
        elif subscripts == 'abba->ab':
            out._array[:] += factor*np.einsum('xab,xba->ab', self._array, self._array2)
        return out

    def slice_to_three(self, subscripts, out=None, factor=1.0, clear=True):
        """Returns a three-index contraction of the four-index object.

           **Arguments:**

           subscripts
                Any of ``abcc->bac``, ``abcc->abc``, ``abcb->abc``, ``abbc->abc``

           **Optional arguments:**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`
        """
        # Error checking
        check_options('subscripts', subscripts, 'abcc->bac', 'abcc->abc', 'abcb->abc', 'abbc->abc')
        if out is None:
            out = DenseThreeIndex(self.nbasis)
        else:
            check_type('out', out, DenseThreeIndex)
            if clear:
                out.clear()
        # Actual computation
        if subscripts == 'abbc->abc':
            slice_to_three_abbc_abc(self._array, self._array2, out._array, factor, clear)
        elif subscripts == 'abcc->bac':
            slice_to_three_abcc_bac(self._array, self._array2, out._array, factor, clear)
        elif subscripts == 'abcc->abc':
            slice_to_three_abcc_abc(self._array, self._array2, out._array, factor, clear)
        elif subscripts == 'abcb->abc':
            L_r = np.diagonal(self._array2, axis1=1, axis2=2)
            out._array[:] += factor*np.tensordot(self._array, L_r, [(0,),(0,)]).swapaxes(1,2)
        return out

    def contract_two_to_four(self, subscripts, two, out=None, factor=1.0, clear=True):
        '''Contracts with a two-index object to obtain a four-index object.

           **Arguments:**

           subscripts
                Any of ``abcd,cd->acbd``, ``abcd,cd->acdb``, ``abcd,cb->acdb``,
                ``abcd,cb->acbd``, ``abcd,ab->acbd``, ``abcd,ab->acdb``,
                ``abcd,ad->acbd``, ``abcd,ad->acdb``, ``abcd,ad->abcd``,
                ``abcd,ad->abdc``, ``abcd,bd->abcd``, ``abcd,bd->abdc``,
                ``abcd,bc->abdc``, ``abcd,bc->abcd``, ``abcd,ac->abcd``,
                ``abcd,ac->abdc``

           two
                An instance of DenseTwoIndex.

           **Optional arguments:**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`
        '''
        check_options('subscripts', subscripts, 'abcd,cd->acbd',
            'abcd,cd->acdb', 'abcd,cb->acdb', 'abcd,cb->acbd', 'abcd,ab->acbd',
            'abcd,ab->acdb', 'abcd,ad->acbd', 'abcd,ad->acdb', 'abcd,ad->abcd',
            'abcd,ad->abdc', 'abcd,bd->abcd', 'abcd,bd->abdc', 'abcd,bc->abdc',
            'abcd,bc->abcd', 'abcd,ac->abcd', 'abcd,ac->abdc')
        raise NotImplementedError

    def contract_two_to_two(self, subscripts, two, out=None, factor=1.0, clear=True):
        """Contract self with a two-index to obtain a two-index.

           **Arguments:**

           subscripts
                Any of ``abcd,bd->ac`` (direct), ``abcd,cb->ad`` (exchange)

           two
                The input two-index object. (DenseTwoIndex)

           **Optional arguments:**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`
        """
        check_options('subscripts', subscripts, 'abcd,bd->ac', 'abcd,cb->ad')
        if out is None:
            out = DenseTwoIndex(self.nbasis)
            if clear:
                out.clear()
        else:
            check_type('out', out, DenseTwoIndex)
        if subscripts == 'abcd,bd->ac':
            tmp = np.tensordot(self._array2, two._array, axes=([(1,2),(1,0)]))
            out._array[:] += factor*np.tensordot(self._array, tmp, [0,0])
        elif subscripts == 'abcd,cb->ad':
            tmp = np.tensordot(self._array2, two._array, axes=([1,1]))
            out._array[:] += factor*np.tensordot(self._array, tmp, ([0,2],[0,2]))
        return out

    def assign_four_index_transform(self, ao_integrals, exp0, exp1=None, exp2=None, exp3=None, method='tensordot'):
        '''Perform four index transformation.

           **Arguments:**

           oa_integrals
                A CholeskyFourIndex with integrals in atomic orbitals.

           exp0
                A DenseExpansion object with molecular orbitals

           **Optional arguments:**

           exp1, exp2, exp3
                Can be provided to transform each index differently. See
                ``parse_four_index_transform_exps`` for details.

           method
                Either ``einsum`` or ``tensordot`` (default).
        '''
        check_type('ao_integrals', ao_integrals, CholeskyFourIndex)
        exp0, exp1, exp2, exp3 = parse_four_index_transform_exps(exp0, exp1, exp2, exp3, DenseExpansion)
        if method == 'einsum':
            if ao_integrals.is_decoupled or not (exp0 is exp1 and exp2 is exp3):
                self.decouple_array2()
                self._array2[:] = np.einsum('bi,kbd->kid', exp1.coeffs, ao_integrals._array2)
                self._array2[:] = np.einsum('dj,kid->kij', exp3.coeffs, self._array2)
            self._array[:] = np.einsum('ai,kac->kic', exp0.coeffs, ao_integrals._array)
            self._array[:] = np.einsum('cj,kic->kij', exp2.coeffs, self._array)
        elif method == 'tensordot':
            if ao_integrals.is_decoupled or not (exp0 is exp1 and exp2 is exp3):
                self.decouple_array2()
                self._array2[:] = np.tensordot(ao_integrals._array2, exp1.coeffs, axes=([1],[0]))
                self._array2[:] = np.tensordot(self._array2, exp3.coeffs, axes=([1],[0]))
            self._array[:] = np.tensordot(ao_integrals._array, exp0.coeffs, axes=([1],[0]))
            self._array[:] = np.tensordot(self._array, exp2.coeffs, axes=([1],[0]))
        else:
            raise ValueError('The method must either be \'einsum\' or \'tensordot\'.')
