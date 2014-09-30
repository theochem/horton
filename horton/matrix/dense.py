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
"""Dense matrix implementations


   Naming scheme for the contract, slice or expand methods
   -------------------------------------------------------

   The name of ``contract``, ``slice`` and ``expand`` methods is as follows::

        [iadd_]{contract,slice,expand}[_X][_Y][_to_Z]

   where each part between square brackets is optional. ``X``, ``Y`` and ``Z``
   can be any of ``one``, ``two``, ``three`` or ``four``. The names
   ``contract``, ``slice`` and ``expand`` have the following meaning:

   ``contract``
        Some actual contraction takes place, i.e. summing of (products) of
        matrix elements.

   ``slice``
        A subset of the elements is merely selected, without any addition of
        elements. This is never combined with ``iadd_``.

   ``expand``
        Products of elements are computed but these products are not added.
        Similar to an outer product but more general.

   When ``iadd_`` is used as a prefix, the result of the contraction is added
   in-place to the self object. In this case, the ``_to_Z`` part is never
   present. A contraction of all input arguments is made. The dimensionality
   of the input arguments is indicated by ``_X`` and ``_Y``.

   When ``_to_Z`` is present, the contraction involves self and possibly other
   arguments whose dimensionality is indicated by ``_X`` and ``_Y``. In this
   case, ``iadd_`` can not be present. The result is written to an output
   argument. If the output argument is not provided, fresh memory is allocated
   to compute the contraction and the result is returned. (This allocation is
   provided for convenience but should not be used in critical situations.)

   Some remarks:

   * Similar conventions apply to an ``expand`` method.
   * All ``contract`` and ``expand`` methods are implemented with the driver
     method ``DenseLinalgFactory.einsum``. However, other implementations than
     `Dense` are free to implement things differently.
   * All ``contract`` and ``expand`` methods never touch the internals of
     higher-index objects.

   For more details, see the docstrings.
"""


import numpy as np
from scipy.linalg import eigh
from horton.log import log
from horton.utils import check_type, check_options
from horton.matrix.base import LinalgFactory, LinalgObject, OneIndex, \
    Expansion, TwoIndex, ThreeIndex, FourIndex, parse_four_index_transform_exps


__all__ = [
    'DenseLinalgFactory', 'DenseOneIndex', 'DenseExpansion',
    'DenseTwoIndex', 'DenseThreeIndex', 'DenseFourIndex'
]


class DenseLinalgFactory(LinalgFactory):
    #
    # DenseOneIndex constructor with default arguments
    #

    def create_one_index(self, nbasis=None):
        '''Create a DenseOneIndex with defaults from the LinalgFactory

           **Optional arguments:**

           nbasis
                The number of basis functions. When not given, the
                default_nbasis value of the DenseLinalgFactory instance will be
                used.
        '''
        nbasis = nbasis or self.default_nbasis
        return DenseOneIndex(nbasis)

    def _check_one_index_init_args(self, other, nbasis=None):
        '''Is an object is compatible the constructor arguments?'''
        nbasis = nbasis or self.default_nbasis
        other.__check_init_args__(nbasis)

    create_one_index.__check_init_args__ = _check_one_index_init_args

    #
    # DenseExpansion constructor with default arguments
    #

    def create_expansion(self, nbasis=None, nfn=None):
        '''Create a DenseExpansion with defaults from the LinalgFactory

           **Optional arguments:**

           nbasis
                The number of basis functions. When not given, the
                default_nbasis value of the DenseLinalgFactory instance will be
                used.

           nfn
                The number of orbitals. When not given, the default_nbasis
                value of the DenseLinalgFactory instance will be used.
        '''
        nbasis = nbasis or self.default_nbasis
        nfn = nfn or nbasis
        return DenseExpansion(nbasis, nfn)

    def _check_expansion_init_args(self, expansion, nbasis=None, nfn=None):
        '''Is an object is compatible the constructor arguments?'''
        nbasis = nbasis or self.default_nbasis
        nfn = nfn or nbasis
        expansion.__check_init_args__(nbasis, nfn)

    create_expansion.__check_init_args__ = _check_expansion_init_args

    #
    # DenseTwoIndex constructor with default arguments
    #

    def create_two_index(self, nbasis=None, nfn=None):
        '''Create a DenseTwoIndex with defaults from the LinalgFactory

           **Optional arguments:**

           nbasis
                The number of basis functions. When not given, the
                default_nbasis value of the DenseLinalgFactory instance will be
                used.

           nfn
                The number of foo. When not given, the default_nbasis
                value of the DenseLinalgFactory instance will be used.
        '''
        nbasis = nbasis or self.default_nbasis
        # Don't replace nfn by self.default_nbasis when it is None! It is
        # a genuine optional argument.
        return DenseTwoIndex(nbasis, nfn)

    def _check_two_index_init_args(self, other, nbasis=None, nfn=None):
        '''Is an object is compatible the constructor arguments?'''
        nbasis = nbasis or self.default_nbasis
        # Don't replace nfn by self.default_nbasis when it is None! It is
        # a genuine optional argument.
        other.__check_init_args__(nbasis, nfn)

    create_two_index.__check_init_args__ = _check_two_index_init_args

    #
    # DenseThreeIndex constructor with default arguments
    #

    def create_three_index(self, nbasis=None):
        '''Create a DenseThreeIndex with defaults from the LinalgFactory

           **Optional arguments:**

           nbasis
                The number of basis functions. When not given, the
                default_nbasis value of the DenseLinalgFactory instance will be
                used.
        '''
        nbasis = nbasis or self.default_nbasis
        return DenseThreeIndex(nbasis)

    def _check_three_index_init_args(self, other, nbasis=None):
        '''Is an object is compatible the constructor arguments?'''
        nbasis = nbasis or self.default_nbasis
        other.__check_init_args__(nbasis)

    create_three_index.__check_init_args__ = _check_three_index_init_args

    #
    # DenseFourIndex constructor with default arguments
    #

    def create_four_index(self, nbasis=None):
        '''Create a DenseFourIndex with defaults from the LinalgFactory

           **Optional arguments:**

           nbasis
                The number of basis functions. When not given, the
                default_nbasis value of the DenseLinalgFactory instance will be
                used.
        '''
        nbasis = nbasis or self.default_nbasis
        return DenseFourIndex(nbasis)

    def _check_four_index_init_args(self, other, nbasis=None):
        '''Is an object is compatible the constructor arguments?'''
        nbasis = nbasis or self.default_nbasis
        other.__check_init_args__(nbasis)

    create_four_index.__check_init_args__ = _check_four_index_init_args

    #
    # Other code
    #

    @staticmethod
    def _allocate_check_output(out, outshape):
        '''Allocate/check the output for the wrappers below

           **Arguments:**

           out
                The output argument.

           outshape
                The expected shape of the output argument.

           **Returns:** the output argument.
        '''
        if len(outshape) == 3 and not (outshape[0] == outshape[1] and
                                       outshape[0] == outshape[2]):
            raise TypeError('A 3-index object must have the same size in all indexes.')
        if len(outshape) == 4 and not (outshape[0] == outshape[1] and
                                       outshape[0] == outshape[2] and
                                       outshape[0] == outshape[3]):
            raise TypeError('A 4-index object must have the same size in all indexes.')

        # Handle the output argument
        if out is None:
            if len(outshape) == 0:
                pass
            elif len(outshape) == 1:
                out = DenseOneIndex(outshape[0])
            elif len(outshape) == 2:
                out = DenseTwoIndex(*outshape)
            elif len(outshape) == 3:
                out = DenseThreeIndex(outshape[0])
            elif len(outshape) == 4:
                out = DenseFourIndex(outshape[0])
            else:
                raise TypeError('The outshape must have length 0, 1, 2, 3 or 4.')
        else:
            if len(outshape) == 0:
                raise TypeError('No output argument can be given when contracting to a scalar.')
            elif len(outshape) == 1:
                check_type('out', out, DenseOneIndex)
            elif len(outshape) == 2:
                check_type('out', out, DenseTwoIndex)
            elif len(outshape) == 3:
                check_type('out', out, DenseThreeIndex)
            elif len(outshape) == 4:
                check_type('out', out, DenseFourIndex)
            else:
                raise TypeError('The outshape must have length 0, 1, 2, 3 or 4.')
            if out.shape != outshape:
                raise TypeError('The output argument does not have the right shape')
        return out

    @staticmethod
    def einsum(subscripts, out, factor, clear, *operands):
        '''einsum wrapper for DenseNBody objects

           **Arguments:**

           subscripts
                The usual subscripts argument for the einsum function, except
                that the format is slightly limited. All indexes must be
                explicitly stated. If the output argument is a scalar, ``->``
                may not be present. In all other cases, ``->`` must be present.

           out
                An output DenseNBody object. This argument is not allowed when
                the result is a scalar.

           factor
                Multiply the contraction by a scalar factor.

           clear
                When set to False, the output argument is not zeroed before
                adding the contraction.

           operands
                A list where each element is (i) a DenseNBody object or (ii) a
                two-tuple of DenseNBody object and ranges to select for slicing.
                The ranges must be a tuple of integers of the form (begin0,
                end0, begin1, end1, ...) where the number of (begin, end) pairs
                depends on the dimensionality of the corresponding tensor.

           **Returns:** the out object. When the out argument is given, this is
           returned with in-place modifications.
        '''
        check_type('subscripts', subscripts, str)
        check_type('factor', factor, float, int)
        check_type('clear', clear, bool)

        # Basic checking of operands and adding ranges if omitted
        new_operands = []
        for operand in operands:
            if isinstance(operand, tuple):
                if len(operand) != 2:
                    raise TypeError('When an operand is a tuple, it must count two elements.')
                tensor, ranges = operand
                if len(tensor.shape)*2 != len(ranges):
                    raise TypeError('The dimensionality of the tensor and the size of ranges does not match.')
            else:
                tensor = operand
                ranges = sum([(0, s) for s in tensor.shape], ())
            new_operands.append((tensor, list(ranges)))
        operands = new_operands

        # Impose explicit format for subscripts and do some basic checking with
        # operands
        if subscripts.count('->') == 1:
            inscripts, outscript = subscripts.split('->')
        elif subscripts.count('->') == 0:
            inscripts = subscripts
            outscript = ''
        else:
            raise ValueError('The subscripts argument must contain a single or none \'->\'.')
        inscripts = inscripts.split(',')
        for inscript, (tensor, ranges) in zip(inscripts, operands):
            if len(inscript) == 1:
                if not isinstance(tensor, DenseOneIndex):
                    raise TypeError('Expecting a DenseOneIndex tensor for subscripts \'%s\'' % inscript)
            elif len(inscript) == 2:
                if not isinstance(tensor, DenseTwoIndex):
                    raise TypeError('Expecting a DenseTwoIndex tensor for subscripts \'%s\'' % inscript)
            elif len(inscript) == 3:
                if not isinstance(tensor, DenseThreeIndex):
                    raise TypeError('Expecting a DenseThreeIndex tensor for subscripts \'%s\'' % inscript)
            elif len(inscript) == 4:
                if not isinstance(tensor, DenseFourIndex):
                    raise TypeError('Expecting a DenseFourIndex tensor for subscripts \'%s\'' % inscript)
            else:
                raise ValueError('The number of subscripts for one tensor must be 1, 2, 3 or 4')

        # take slices of the input arrays
        new_operands = []
        for tensor, ranges in operands:
            slices = []
            while len(ranges) > 0:
                begin = ranges.pop(0)
                end = ranges.pop(0)
                slices.append(slice(begin, end))
            array = tensor._array[tuple(slices)]
            new_operands.append(array)
        operands = new_operands

        # Determine shape of the output
        outshape = []
        for outchar in outscript:
            size = None
            for inscript, array in zip(inscripts, operands):
                for i, inchar in enumerate(inscript):
                    if inchar == outchar:
                        if size is None:
                            size = array.shape[i]
                        else:
                            if size != array.shape[i]:
                                raise TypeError('Incompatible shapes of input operands.')
            if size is None:
                raise TypeError('Could not determine the shape of the output.')
            outshape.append(size)
        outshape = tuple(outshape)

        # Allocate/check the output argument
        out = DenseLinalgFactory._allocate_check_output(out, outshape)

        # do the actual work
        if len(outshape) == 0:
            out = float(np.einsum(subscripts + '->...', *operands))*factor
            assert isinstance(out, float)
        else:
            if clear:
                # can be done without temporary
                out._array[:] = np.einsum(subscripts, *operands)
                # FIXME: due to a bug in numpy, we can't use the output argument
                # for all cases. For more details, see:
                # https://github.com/numpy/numpy/issues/5147
                #np.einsum(subscripts, *operands, out=out._array)
                out.iscale(factor)
            else:
                # can't be done without temporary
                out._array += np.einsum(subscripts, *operands)*factor
        return out

    @staticmethod
    def tensordot(a, b, axes, out=None, factor=1.0, clear=True):
        '''tensordot wrapper for dense n-index objects.

           **Arguments:**

           a, b, axes
                See documentation of numpy.tensordot

           **Optional arguments:**

           out
                The output argument, an instance of one of the Dense?Index classes.
        '''
        # Basic type checking
        check_type('a', a, DenseOneIndex, DenseTwoIndex, DenseThreeIndex, DenseFourIndex)
        check_type('b', b, DenseOneIndex, DenseTwoIndex, DenseThreeIndex, DenseFourIndex)

        # Determine the shape of the output
        outshape = tuple([a.shape[i] for i in xrange(a.ndim) if i not in axes[0]]) + \
                   tuple([b.shape[i] for i in xrange(b.ndim) if i not in axes[1]])

        # Allocate/check output
        out = DenseLinalgFactory._allocate_check_output(out, outshape)

        # Actual work
        if len(outshape) == 0:
            out = np.tensordot(a._array, b._array, axes)*factor
        else:
            if clear:
                # can't be done without temporary
                out._array[:] = np.tensordot(a._array, b._array, axes)
                out.iscale(factor)
            else:
                # can't be done without temporary
                out._array += np.tensordot(a._array, b._array, axes)*factor
        return out


class DenseOneIndex(OneIndex):
    """Dense one-dimensional matrix (vector)

       This is also used for (diagonal) density matrices.
    """

    #
    # Constructor and destructor
    #

    def __init__(self, nbasis):
        """
           **Arguments:**

           nbasis
                The number of basis functions.
        """
        self._array = np.zeros((nbasis,))
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)

    #
    # Methods from base class
    #

    def __check_init_args__(self, nbasis):
        '''Is self compatible with the given constructor arguments?'''
        assert nbasis == self.nbasis

    def __eq__(self, other):
        '''Compare self with other'''
        return isinstance(other, DenseOneIndex) and \
            other.nbasis == self.nbasis and \
            (other._array == self._array).all()

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct an instance from data previously stored in an h5py.Group.

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        nbasis = grp['array'].shape[0]
        result = cls(nbasis)
        grp['array'].read_direct(result._array)
        return result

    def to_hdf5(self, grp):
        '''Dump this object in an h5py.Group

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        grp.attrs['class'] = self.__class__.__name__
        grp['array'] = self._array

    def new(self):
        '''Return a new one-index object with the same nbasis'''
        return DenseOneIndex(self.nbasis)

    def _check_new_init_args(self, other):
        '''Check whether an already initialized object is compatible'''
        other.__check_init_args__(self.nbasis)

    new.__check_init_args__ = _check_new_init_args

    def clear(self):
        '''Reset all elements to zero.'''
        self._array[:] = 0.0

    def copy(self, begin=0, end=None):
        '''Return a copy of (a part of) the object

           **Optional arguments:**

           begin, end
                Can be used to select a subblock of the object. When not given,
                the full range is used.
        '''
        if end is None:
            end = self.nbasis
        result = DenseOneIndex(end - begin)
        result._array[:] = self._array[begin:end]
        return result

    def assign(self, other):
        '''Assign with the contents of another object

           **Arguments:**

           other
                Another DenseOneIndex object or a scalar.
        '''
        if isinstance(other, DenseOneIndex):
            self._array[:] = other._array
        elif isinstance(other, float) or isinstance(other, int):
            self._array[:] = other
        else:
            raise TypeError('Do not know how to assign object of type %s.' % type(other))

    def randomize(self):
        '''Fill with random normal data'''
        self._array[:] = np.random.normal(0, 1, self.shape)

    def permute_basis(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.

           **Arguments:**

           permutation
                An integer numpy array that defines the new order of the basis
                functions.
        '''
        self._array[:] = self._array[permutation]

    def change_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.

           **Arguments:**

           signs
                A numpy array with sign changes indicated by +1 and -1.
        '''
        self._array *= signs

    def iadd(self, other, factor=1.0):
        '''Add another DenseOneIndex object in-place, multiplied by factor

           **Arguments:**

           other
                A DenseOneIndex instance to be added

           **Optional arguments:**

           factor
                A scalar factor
        '''
        check_type('other', other, DenseOneIndex)
        check_type('factor', factor, float, int)
        self._array += other._array*factor

    def iscale(self, factor):
        '''In-place multiplication with a scalar

           **Arguments:**

           factor
                A scalar factor.
        '''
        check_type('factor', factor, float, int)
        self._array *= factor

    def get_element(self, i):
        '''Return a matrix element'''
        return self._array[i]

    def set_element(self, i, value):
        '''Set a matrix element'''
        self._array[i] = value

    #
    # Properties
    #

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._array.shape[0]

    nbasis = property(_get_nbasis)

    def _get_shape(self):
        '''The shape of the object'''
        return self._array.shape

    shape = property(_get_shape)


class DenseExpansion(Expansion):
    """An expansion of several functions in a basis with a dense matrix of
       coefficients. The implementation is such that the columns of self._array
       contain the orbitals.
    """

    #
    # Constructor and destructor
    #

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

    #
    # Methods from base class
    #

    def __check_init_args__(self, nbasis, nfn=None):
        '''Is self compatible with the given constructor arguments?'''
        if nfn is None:
            nfn = nbasis
        assert nbasis == self.nbasis
        assert nfn == self.nfn

    def __eq__(self, other):
        '''Compare self with other'''
        return isinstance(other, DenseExpansion) and \
            other.nbasis == self.nbasis and \
            other.nfn == self.nfn and \
            (other._coeffs == self._coeffs).all() and \
            (other._energies == self._energies).all() and \
            (other._occupations == self._occupations).all()

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct an instance from data previously stored in an h5py.Group.

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        if grp.attrs['class'] != cls.__name__:
            raise TypeError('The class of the expansion in the HDF5 file does not match.')
        nbasis, nfn = grp['coeffs'].shape
        result = cls(nbasis, nfn)
        grp['coeffs'].read_direct(result._coeffs)
        grp['energies'].read_direct(result._energies)
        grp['occupations'].read_direct(result._occupations)
        return result

    def to_hdf5(self, grp):
        '''Dump this object in an h5py.Group

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        grp.attrs['class'] = self.__class__.__name__
        grp['coeffs'] = self._coeffs
        grp['energies'] = self._energies
        grp['occupations'] = self._occupations

    def new(self):
        '''Return a new expansion object with the same nbasis and nfn'''
        return DenseExpansion(self.nbasis, self.nfn)

    def _check_new_init_args(self, other):
        '''Check whether an already initialized object is compatible'''
        other.__check_init_args__(self.nbasis, self.nfn)

    new.__check_init_args__ = _check_new_init_args

    def clear(self):
        '''Reset all elements to zero.'''
        self._coeffs[:] = 0.0
        self._energies[:] = 0.0
        self._occupations[:] = 0.0

    def copy(self):
        '''Return a copy of the object'''
        result = DenseExpansion(self.nbasis, self.nfn)
        result._coeffs[:] = self._coeffs
        result._energies[:] = self._energies
        result._occupations[:] = self._occupations
        return result

    def assign(self, other):
        '''Assign with the contents of another object

           **Arguments:**

           other
                Another DenseExpansion object.
        '''
        check_type('other', other, DenseExpansion)
        self._coeffs[:] = other._coeffs
        self._energies[:] = other._energies
        self._occupations[:] = other._occupations

    def randomize(self):
        '''Fill with random normal data'''
        self._coeffs[:] = np.random.normal(0, 1, self._coeffs.shape)
        self._energies[:] = np.random.normal(0, 1, self._energies.shape)
        self._occupations[:] = np.random.normal(0, 1, self._occupations.shape)

    def permute_basis(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.

           **Arguments:**

           permutation
                An integer numpy array that defines the new order of the basis
                functions.
        '''
        self._coeffs[:] = self.coeffs[permutation]

    def change_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.

           **Arguments:**

           signs
                A numpy array with sign changes indicated by +1 and -1.
        '''
        self._coeffs *= signs.reshape(-1,1)

    def check_normalization(self, overlap, eps=1e-4):
        '''Run an internal test to see if the orbitals are normalized

           **Arguments:**

           overlap
                The overlap two-index operator

           **Optional arguments:**

           eps
                The allowed deviation from unity, very loose by default.
        '''
        check_type('overlap', overlap, DenseTwoIndex)
        for i in xrange(self.nfn):
            if self.occupations[i] == 0:
                continue
            norm = overlap.inner(self._coeffs[:,i], self._coeffs[:,i])
            #print i, norm
            assert abs(norm-1) < eps, 'The orbitals are not normalized!'

    def error_eigen(self, fock, overlap):
        """Compute the error of the orbitals with respect to the eigenproblem

           **Arguments:**

           fock
                A DenseTwoIndex Hamiltonian (or Fock) operator.

           overlap
                A DenseTwoIndex overlap operator.

           **Returns:** the RMSD error on the orbital energies
        """
        check_type('fock', fock, DenseTwoIndex)
        check_type('overlap', overlap, DenseTwoIndex)
        errors = np.dot(fock._array, (self.coeffs)) \
                 - self.energies*np.dot(overlap._array, (self.coeffs))
        return np.sqrt((abs(errors)**2).mean())

    #
    # Properties
    #

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
        return self._energies.view()

    energies = property(_get_energies)

    def _get_occupations(self):
        '''The orbital occupations'''
        return self._occupations.view()

    occupations = property(_get_occupations)

    #
    # New methods for this implementation
    # TODO: consider adding these to base class
    #

    def from_fock(self, fock, overlap):
        '''Diagonalize a Fock matrix to obtain orbitals and energies

           **Arguments:**

           fock
                The fock matrix, an instance of DenseTwoIndex.

           overlap
                The overlap matrix, an instance of DenseTwoIndex.
        '''
        check_type('fock', fock, DenseTwoIndex)
        check_type('overlap', overlap, DenseTwoIndex)
        evals, evecs = eigh(fock._array, overlap._array)
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
        check_type('dm', dm, DenseTwoIndex)
        check_type('overlap', overlap, DenseTwoIndex)
        # Construct a level-shifted operator
        occ = overlap.copy()
        occ.idot(dm)
        occ.idot(overlap)
        # diagonalize and compute eigenvalues
        evals, evecs = eigh(occ._array, overlap._array)
        self._coeffs[:] = evecs[:,:self.nfn]
        self._occupations[:] = evals
        self._energies[:] = 0.0

    def get_homo_index(self, offset=0):
        '''Return the index of a HOMO orbital.

           **Optional arguments**:

           offset
                By default, the (highest) homo energy is returned. When this
                index is above zero, the corresponding lower homo energy is
                returned.
        '''
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
        '''Return the index of a LUMO orbital.

           **Optional arguments**:

           offset
                By default, the (lowest) lumo energy is returned. When this
                index is above zero, the corresponding higher homo energy is
                returned.
        '''
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

    def to_dm(self, out=None, factor=1.0, clear=True, other=None):
        """Compute the density matrix

           **Optional arguments:**

           out
                An output density matrix (DenseTwoIndex instance).

           factor
                The density matrix is multiplied by the given scalar.

           clear
                When set to False, the output density matrix is not zeroed first.

           other
                Another DenseExpansion object to construct a transfer-density
                matrix. When combined with an output argument, it may not be
                a Hermitian two-index object.
        """
        # parse first argument
        if out is None:
            if other is None:
                out = DenseTwoIndex(self.nbasis)
            else:
                out = DenseTwoIndex(self.nbasis, self.nbasis)
        else:
            check_type('out', out, DenseTwoIndex)
            if other is not None and out.is_hermitian:
                raise TypeError('A transfer density matrix can only be stored in a non-Hermitian DenseTwoIndex.')

        # actual computation
        if clear:
            if other is None:
                out._array[:] = factor*np.dot(self._coeffs*self.occupations, self._coeffs.T)
            else:
                out._array[:] = factor*np.dot(self._coeffs*self.occupations, other._coeffs.T)
        else:
            if other is None:
                out._array += factor*np.dot(self._coeffs*self.occupations, self._coeffs.T)
            else:
                out._array += factor*np.dot(self._coeffs*self.occupations, other._coeffs.T)
        return out

    def assign_dot(self, other, tf2):
        '''Dot product of orbitals in a DenseExpansion and TwoIndex object

           **Arguments:**

           other
                An expansion object with input orbitals

           tf2
                A two-index object

           The transformed orbitals are stored in self.
        '''
        check_type('other', other, DenseExpansion)
        check_type('tf2', tf2, DenseTwoIndex)
        if not (self.nbasis == other.nbasis):
            raise TypeError('Both expansions must have the same number of basis functions.')
        if not (tf2.shape[0] == other.nfn and tf2.shape[1] == self.nfn):
            raise TypeError('The shape of the two-index object is incompatible with that of the expansions.')
        self._coeffs[:] = np.dot(other.coeffs, tf2._array)

    def rotate_random(self):
        '''Apply random unitary transformation distributed with Haar measure

           The attributes ``energies`` and ``occupations`` are not altered.
        '''
        z = np.random.normal(0, 1, (self.nfn, self.nfn))
        q, r = np.linalg.qr(z)
        self.coeffs[:] = np.dot(self.coeffs, q)

    def rotate_2orbitals(self, angle=0.7853981633974483, index0=None, index1=None):
        '''Rotate two orbitals

           **Optional arguments:**

           angle
                The rotation angle, defaults to 45 deg.

           index0, index1
                The orbitals to rotate, defaults to HOMO and LUMO,

           The attributes ``energies`` and ``occupations`` are not altered.
        '''
        if index0 == None:
            index0 = self.get_homo_index()
        if index1 == None:
            index1 = self.get_lumo_index()
        old0 = self.coeffs[:,index0].copy()
        old1 = self.coeffs[:,index1].copy()
        self.coeffs[:,index0] = np.cos(angle)*old0 - np.sin(angle)*old1
        self.coeffs[:,index1] = np.sin(angle)*old0 + np.cos(angle)*old1

    def swap_orbitals(self, swaps):
        '''Change the order of the orbitals using pair-exchange

           **Arguments:**

           swaps
                An integer numpy array with two columns where every row
                corresponds to one swap.

           The attributes ``energies`` and ``occupations`` are also reordered.
        '''
        if not (swaps.shape[1] == 2 and swaps.ndim == 2 and issubclass(swaps.dtype.type, int)):
            raise TypeError('The argument swaps has the wrong shape/type.')
        for iswap in range(len(swaps)):
            index0, index1 = swaps[iswap]
            if log.do_medium:
                log('  Swapping orbitals %i and %i' %(index0, index1))
            tmp = self.coeffs[:,index0].copy()
            self.coeffs[:,index0] = self.coeffs[:,index1]
            self.coeffs[:,index1] = tmp
            self.energies[index0], self.energies[index1] =\
                self.energies[index1], self.energies[index0]
            self.occupations[index0], self.occupations[index1] =\
                self.occupations[index1], self.occupations[index0]


class DenseTwoIndex(TwoIndex):
    """Dense symmetric two-dimensional matrix, also used for density matrices.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """

    #
    # Constructor and destructor
    #

    def __init__(self, nbasis, nfn=None):
        """
           **Arguments:**

           nbasis
                The number of basis functions. (Number of rows. Also number of
                columns, unless nfn is given.)

           **Optional arguments:**

           nfn
                When given, this is the number of columns (second index).

           Note that by default the two-index object is assumed to be Hermitian.
           Only when nfn is given, this assumption is dropped.
        """
        if nfn is None:
            nfn = nbasis
            self._is_hermitian = True
        else:
            self._is_hermitian = False
        self._array = np.zeros((nbasis, nfn))
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)

    #
    # Methods from base class
    #

    def __check_init_args__(self, nbasis, nfn=None):
        '''Is self compatible with the given constructor arguments?'''
        assert nbasis == self.nbasis
        assert nfn == self.nfn

    def __eq__(self, other):
        '''Compare self with other'''
        return isinstance(other, DenseTwoIndex) and \
            other.nbasis == self.nbasis and \
            other.nfn == self.nfn and \
            (other._array == self._array).all()

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct an instance from data previously stored in an h5py.Group.

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        nbasis = grp['array'].shape[0]
        if grp.attrs['is_hermitian']:
            nfn = None
        else:
            nfn = grp['array'].shape[1]
        result = cls(nbasis, nfn)
        grp['array'].read_direct(result._array)
        return result

    def to_hdf5(self, grp):
        '''Dump this object in an h5py.Group

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        grp.attrs['class'] = self.__class__.__name__
        grp.attrs['is_hermitian'] = self._is_hermitian
        grp['array'] = self._array

    def new(self):
        '''Return a new two-index object with the same nbasis (and nfn)'''
        return DenseTwoIndex(self.nbasis, self.nfn)

    def _check_new_init_args(self, other):
        '''Check whether an already initialized object is compatible'''
        other.__check_init_args__(self.nbasis, self.nfn)

    new.__check_init_args__ = _check_new_init_args

    def clear(self):
        '''Reset all elements to zero.'''
        self._array[:] = 0.0

    def copy(self, begin0=0, end0=None, begin1=0, end1=None):
        '''Return a copy of (a part of) the object

           **Optional arguments:**

           begin0, end0, begin1, end1
                Can be used to select a subblock of the object. When not given,
                the full range is used. Only when the ranges are equal for both
                axes, an Hermitian two-index may be returned.
        '''
        end0, end1 = self._fix_ends(end0, end1)
        nbasis = end0 - begin0
        if self._is_hermitian and begin0 == begin1 and end0 == end1:
            nfn = None
        else:
            nfn = end1 - begin1
        result = DenseTwoIndex(nbasis, nfn)
        result._array[:] = self._array[begin0:end0, begin1:end1]
        return result

    def assign(self, other):
        '''Assign a new contents to the two-index object

           **Arguments:**

           other
                The new data, may be DenseTwoIndex or a scalar value.
        '''
        if isinstance(other, DenseTwoIndex):
            self._array[:] = other._array
        elif isinstance(other, float) or isinstance(other, int):
            self._array[:] = other
        else:
            raise TypeError('Do not know how to assign object of type %s.' % type(other))

    def iadd(self, other, factor=1.0, begin0=0, end0=None, begin1=0, end1=None, transpose=False):
        '''Add another DenseTwoIndex object in-place, multiplied by factor

           **Arguments:**

           other
                A DenseTwoIndex instance to be added

           **Optional arguments:**

           factor
                A scalar factor

           begin0, end0, begin1, end1
                When given, specify the ranges where the contribution will be
                added. When not given, the full range is used.
        '''
        check_type('factor', factor, float, int)
        end0, end1 = self._fix_ends(end0, end1)
        if isinstance(other, DenseTwoIndex):
            if transpose:
                self._array[begin0:end0, begin1:end1] += other._array.T*factor
            else:
                self._array[begin0:end0, begin1:end1] += other._array*factor
        elif isinstance(other, float):
            self._array[begin0:end0, begin1:end1] += other*factor
        else:
            raise TypeError('Do not know how to add in-place an object of type %s.' % type(other))

    def iscale(self, factor):
        '''In-place multiplication with a scalar

           **Arguments:**

           factor
                A scalar factor.
        '''
        check_type('factor', factor, float, int)
        self._array *= factor

    def randomize(self):
        '''Fill with random normal data'''
        self._array[:] = np.random.normal(0, 1, self.shape)
        if self._is_hermitian:
            self._array += self._array.T
            self._array *= 0.5

    def permute_basis(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.

           **Arguments:**

           permutation
                An integer numpy array that defines the new order of the basis
                functions.
        '''
        self._array[:] = self._array[permutation]
        if self._is_hermitian:
            self._array[:] = self._array[:, permutation]

    def change_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.

           **Arguments:**

           signs
                A numpy array with sign changes indicated by +1 and -1.
        '''
        self._array *= signs
        if self._is_hermitian:
            self._array *= signs.reshape(-1, 1)

    def get_element(self, i, j):
        '''Return a matrix element'''
        return self._array[i, j]

    def set_element(self, i, j, value):
        '''Set a matrix element

           If the two-index object is Hermitian, also the (j,i) element will
           be set appropriately.
        '''
        self._array[i, j] = value
        if self._is_hermitian:
            self._array[j, i] = value

    def sum(self, begin0=0, end0=None, begin1=0, end1=None):
        '''Return the sum of all elements (in the selected range)

           **Optional arguments:**

           begin0, end0, begin1, end1
                Can be used to select a subblock of the object to be contracted.
        '''
        end0, end1 = self._fix_ends(end0, end1)
        return self._array[begin0:end0, begin1:end1].sum()

    def trace(self, begin0=0, end0=None, begin1=0, end1=None):
        '''Return the trace of the two-index object.

           **Optional arguments:**

           begin0, end0, begin1, end1
                Can be used to select a subblock of the object to be contracted.
        '''
        end0, end1 = self._fix_ends(end0, end1)
        if end0-begin0 != end1-begin1:
            raise ValueError('Only the trace of a square (part of a) two-index object can be computed.')
        return np.trace(self._array[begin0:end0, begin1:end1])

    def itranspose(self):
        '''In-place transpose'''
        self._array[:] = self._array.T

    def inner(self, vec0, vec1):
        '''Compute an inner product of two vectors using the two-index as a metric

           **Arguments:**

           vec0, vec1
                The vectors, either DenseOneIndex or numpy arrays.
        '''
        if vec0.shape != (self.shape[0],):
            raise TypeError('The length of vec0 does not match the shape of the two-index object.')
        if vec1.shape != (self.shape[1],):
            raise TypeError('The length of vec1 does not match the shape of the two-index object.')
        if isinstance(vec0, DenseOneIndex) and isinstance(vec1, DenseOneIndex):
            return np.dot(vec0._array, np.dot(self._array, vec1._array))
        elif isinstance(vec0, np.ndarray) and isinstance(vec1, np.ndarray):
            return np.dot(vec0, np.dot(self._array, vec1))
        else:
            raise TypeError('Do not know how to compute inner product with objects of type \'%s\' and \'%s\'.' % (
                type(vec0), type(vec1)
            ))

    #
    # Properties
    #

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._array.shape[0]

    nbasis = property(_get_nbasis)

    def _get_nfn(self):
        '''The other size of the two-index object'''
        if self._is_hermitian:
            return None
        else:
            return self.shape[1]

    nfn = property(_get_nfn)

    def _get_is_hermitian(self):
        '''Whether the two-index object is _supposed_ to be Hermitian.'''
        return self._is_hermitian

    is_hermitian = property(_get_is_hermitian)

    def _get_shape(self):
        '''The shape of the object'''
        return self._array.shape

    shape = property(_get_shape)

    #
    # New methods for this implementation
    # TODO: consider adding these to base class
    #

    def assign_diagonal(self, value):
        '''Set diagonal elements to value

           **Arguments:**

           value
                Either a scalar or a DenseOneIndex object
        '''
        if isinstance(value, DenseOneIndex):
            np.fill_diagonal(self._array, value._array)
        elif isinstance(value, float):
            np.fill_diagonal(self._array, value)
        else:
            raise TypeError('Do not know how to set diagonal with object of type %s.' % type(value))

    def copy_diagonal(self, out=None, begin=0, end=None):
        '''Copy (part of) the diagonal of the two-index object

           **Optional arguments:**

           out
                The output argument (DenseOneIndex with proper size).

           begin, end
                Can be used to select a range of the diagonal. If not given,
                then the entire diagonal is copied.
        '''
        if not self.shape[0] == self.shape[1]:
            raise TypeError('The diagonal can only be copied when the two-index object is squared.')
        if end is None:
            end = self.shape[0]
        if out is None:
            out = DenseOneIndex(end - begin)
        else:
            check_type('out', out, DenseOneIndex)
            if out.shape != (end - begin,):
                raise TypeError('The output argument has the incorrect shape.')
        out._array[:] = np.diagonal(self._array[begin:end,begin:end])
        return out

    def check_symmetry(self):
        '''Check the symmetry of the array. For testing only.'''
        if not self._is_hermitian:
            raise RuntimeError('This makes no sense for a matrix that is not supposed to be Hermitian.')
        assert (self._array == self._array.T).all()

    def symmetrize(self, factor=1.0):
        '''Symmetrize DenseTwoIndex using M_sym=(M+M^\dagger)/2

           **Optional arguments:**

           factor
                Multiply the result by the given factor.
        '''
        check_type('factor', factor, float, int)
        if self.shape[0] != self.shape[1]:
            raise TypeError('Can not symmetrize a rectangular two-index object.')
        result = DenseTwoIndex(self.nbasis)
        result._array[:] = self._array
        result._array += self._array.T
        result.iscale(factor/2.0)
        return result

    def contract_to_one(self, subscripts, out=None, factor=1.0, clear=True, begin0=0, end0=None, begin1=0, end1=None):
        '''Contract self to OneIndex.

           **Arguments:**

           subscripts
                ``ab->b``: contract first index, ``ab->a``: contract second
                index.

           **Optional arguments**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`

           begin0, end0, begin1, end1
                Can be used to contract only a part of the two-index object

           **Returns:** the contracted one-index object.
        '''
        check_options('subscripts', subscripts, 'ab->b', 'ab->a')
        return DenseLinalgFactory.einsum(subscripts, out, factor, clear, (self, (begin0, end0, begin1, end1)))

    def contract_two_to_one(self, subscripts, other, out=None, factor=1.0, clear=True):
        '''Contract two TwoIndex objects to a one OneIndex.

           **Arguments:**

           subscripts:
                ``ab,ab->b``: contract first index. ``ab,ab->a``: contract
                second index.

           other
                The second DenseTwoIndex object.

           **Optional arguments:**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`

           **Returns:** the contracted one-index object.
        '''
        check_options('subscripts', subscripts, 'ab,ab->b', 'ab,ab->a')
        return DenseLinalgFactory.einsum(subscripts, out, factor, clear, self, other)

    def iadd_t(self, other, factor=1.0, begin0=0, end0=None, begin1=0, end1=None):
        '''See :py:meth:`DenseTwoIndex.iadd`, transpose=True'''
        self.iadd(other, factor, begin0, end0, begin1, end1, True)

    def iadd_outer(self, other0, other1, factor=1.0):
        '''In-place addition of outer product of two other DenseTwoIndex

           **Arguments:**

           other0, other1
                Two-index objects that go into the outer product. They are
                raveled before taking the outer product.

           **Optional arguments:**

           factor
                The term added is scaled by this factor.
        '''
        check_type('other0', other0, DenseTwoIndex)
        check_type('other1', other1, DenseTwoIndex)
        check_type('factor', factor, float, int)
        self._array += np.outer(other0._array.ravel(), other1._array.ravel())*factor

    def iadd_kron(self, other0, other1, factor=1.0):
        '''In-place addition of kronecker product of two other DenseTwoIndex

           **Arguments:**

           other0, other1
                Two-index objects that go into the kronecker product.

           **Optional arguments:**

           factor
                The term added is scaled by this factor.
        '''
        check_type('other0', other0, DenseTwoIndex)
        check_type('other1', other1, DenseTwoIndex)
        check_type('factor', factor, float, int)
        self._array += np.kron(other0._array, other1._array)*factor

    def iadd_dot(self, other0, other1, factor=1.0):
        '''In-place addition of dot product: ``other0 * other1``

           **Arguments:**

           other0, other1
                Two-index objects that go into the kronecker product.

           **Optional arguments:**

           factor
                The term added is scaled by this factor.
        '''
        check_type('other0', other0, DenseTwoIndex)
        check_type('other1', other1, DenseTwoIndex)
        check_type('factor', factor, float, int)
        self._array[:] += np.dot(other0._array, other1._array)*factor

    def iadd_tdot(self, other0, other1, factor=1.0):
        '''In-place addition of dot product: ``other0.T * other1``

           **Arguments:**

           other0, other1
                Two-index objects that go into the kronecker product.

           **Optional arguments:**

           factor
                The term added is scaled by this factor.
        '''
        check_type('other0', other0, DenseTwoIndex)
        check_type('other1', other1, DenseTwoIndex)
        check_type('factor', factor, float, int)
        self._array[:] += np.dot(other0._array.T, other1._array)*factor

    def iadd_dott(self, other0, other1, factor=1.0):
        '''In-place addition of dot product: ``other0 * other1.T``

           **Arguments:**

           other0, other1
                Two-index objects that go into the kronecker product.

           **Optional arguments:**

           factor
                The term added is scaled by this factor.
        '''
        check_type('other0', other0, DenseTwoIndex)
        check_type('other1', other1, DenseTwoIndex)
        check_type('factor', factor, float, int)
        self._array[:] += np.dot(other0._array, other1._array.T)*factor

    def iadd_shift(self, lshift):
        '''Add positive shift to elements. If negative replace by shift

           **Arguments:**

           lshift
                A scalar used to augment the matrix elements.
        '''
        check_type('lshift', lshift, float, int)
        self._array[self._array >= 0] += lshift
        self._array[self._array < 0] = lshift

    def iadd_contract_two_one(self, subscripts, two, one, factor=1.0):
        '''In-place addition of a contraction of a two- and one-index object

           **Arguments:**

           two, one
                The two- and one-index objects, respectively. Instances of
                DenseTwoIndex and DenseOneIndex respectively.

           subscripts
                ``ab,b->ab``: contract with the first index of the two-index object.
                ``ab,a->ab``: contract with the second index of the two-index.
                ``ba,a->ab``: contract with the second index and transpose.
                object.

           factor
                The term added is scaled by this factor.
        '''
        check_options('subscripts', subscripts, 'ab,b->ab', 'ab,a->ab', 'ba,a->ab')
        DenseLinalgFactory.einsum(subscripts, self, factor, False, two, one)

    def contract_two(self, subscripts, two):
        '''Compute the trace using with other two-index objects

           **Arguments:**

           subscripts
                ``ab,ba``: trace of matrix product. ``ab,ab``: trace after
                element-wise multiplication (expectation_value).

           two
                A DenseTwoIndex instance
        '''
        check_options('subscripts', subscripts, 'ab,ab', 'ab,ba')
        return DenseLinalgFactory.einsum(subscripts, None, 1.0, True, self, two)

    def idot(self, other):
        '''In-place dot product: self = self * other

           **Arguments:**

           other
                The other array
        '''
        check_type('other', other, DenseTwoIndex)
        self._array[:] = np.dot(self._array, other._array)

    def imul(self, other, factor=1.0):
        '''In-place element-wise multiplication: ``self *= other * factor``

           **Arguments:**

           other
                A DenseTwoIndex instance.

           **Optional arguments:**

           factor
                The two-index object is scaled by this factor.
        '''
        check_type('other', other, DenseTwoIndex)
        check_type('factor', factor, float, int)
        self._array *= other._array
        self.iscale(factor)

    def imul_t(self, other, factor=1.0):
        '''In-place element-wise multiplication: ``self *= other.T * factor``

           **Arguments:**

           other
                A DenseTwoIndex instance.

           **Optional arguments:**

           factor
                The two-index object is scaled by this factor.
        '''
        check_type('other', other, DenseTwoIndex)
        check_type('factor', factor, float, int)
        self._array *= other._array.T
        self.iscale(factor)

    def distance_inf(self, other):
        '''The infinity norm distance between self and other

           **Arguments:**

           other
                A DenseTwoIndex instance.
        '''
        check_type('other', other, DenseTwoIndex)
        return abs(self._array.ravel() - other._array.ravel()).max()

    def iabs(self):
        '''In-place absolute values'''
        self._array[:] = abs(self._array)

    def assign_two_index_transform(self, ao_integrals, exp0, exp1=None):
        '''Perform two index transformation: ``exp0.T * ao_integrals * exp0``

           **Arguments:**

           ao_integrals
                Something computed in the atomic orbital basis

           exp0
                The molecular orbitals.

           **Optional arguments:**

           exp1
                When given, the following is computed: ``exp0.T * ao_integrals *
                exp1``
        '''
        if exp1 is None:
            exp1 = exp0
            self._is_hermitian = ao_integrals.is_hermitian
        else:
            self._is_hermitian = False
        self._array[:] = reduce(np.dot, [exp0.coeffs.T, ao_integrals._array, exp1.coeffs])
        if self._is_hermitian:
            # Matrix product does not respect symmetry.
            # FIXME: Is there a better way to solve this?
            self._array += self._array.T
            self._array *= 0.5


class DenseThreeIndex(ThreeIndex):
    """Dense three-dimensional object.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """

    #
    # Constructor and destructor
    #

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
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)

    #
    # Methods from base class
    #

    def __check_init_args__(self, nbasis):
        '''Is self compatible with the given constructor arguments?'''
        assert nbasis == self.nbasis

    def __eq__(self, other):
        '''Compare self with other'''
        return isinstance(other, DenseThreeIndex) and \
            other.nbasis == self.nbasis and \
            (other._array == self._array).all()

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct an instance from data previously stored in an h5py.Group.

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        nbasis = grp['array'].shape[0]
        result = cls(nbasis)
        grp['array'].read_direct(result._array)
        return result

    def to_hdf5(self, grp):
        '''Dump this object in an h5py.Group

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        grp.attrs['class'] = self.__class__.__name__
        grp['array'] = self._array

    def new(self):
        '''Return a new three-index object with the same nbasis'''
        return DenseThreeIndex(self.nbasis)

    def _check_new_init_args(self, other):
        '''Check whether an already initialized object is compatible'''
        other.__check_init_args__(self.nbasis)

    new.__check_init_args__ = _check_new_init_args

    def clear(self):
        '''Reset all elements to zero.'''
        self._array[:] = 0.0

    def copy(self, begin0=0, end0=None, begin1=0, end1=None, begin2=0, end2=None):
        '''Return a copy of (a part of) the object

           **Optional arguments:**

           begin0, end0, begin1, end1, begin2, end2
                Can be used to select a subblock of the object. When not given,
                the full range is used.
        '''
        end0, end1, end2 = self._fix_ends(end0, end1, end2)
        if end1 - begin1 != end0 - begin0:
            raise TypeError('The shape of the result must be the same along all indexes.')
        if end2 - begin2 != end0 - begin0:
            raise TypeError('The shape of the result must be the same along all indexes.')
        result = DenseThreeIndex(end0 - begin0)
        result._array[:] = self._array[begin0:end0, begin1:end1, begin2:end2]
        return result

    def assign(self, other):
        '''Assign with the contents of another object

           **Arguments:**

           other
                Another DenseThreeIndex object or a scalar.
        '''
        check_type('other', other, DenseThreeIndex)
        self._array[:] = other._array

    def randomize(self):
        '''Fill with random normal data'''
        self._array[:] = np.random.normal(0, 1, self.shape)

    def permute_basis(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.

           **Arguments:**

           permutation
                An integer numpy array that defines the new order of the basis
                functions.
        '''
        self._array[:] = self._array[permutation]
        self._array[:] = self._array[:, permutation]
        self._array[:] = self._array[:, :, permutation]

    def change_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.

           **Arguments:**

           signs
                A numpy array with sign changes indicated by +1 and -1.
        '''
        self._array *= signs
        self._array *= signs.reshape(-1, 1)
        self._array *= signs.reshape(-1, 1, 1)

    def iadd(self, other, factor=1.0):
        '''Add another DenseThreeIndex object in-place, multiplied by factor

           **Arguments:**

           other
                A DenseThreeIndex instance to be added

           **Optional arguments:**

           factor
                The added term is scaled by this factor.
        '''
        check_type('other', other, DenseThreeIndex)
        check_type('factor', factor, float, int)
        self._array += other._array*factor

    def iscale(self, factor):
        '''In-place multiplication with a scalar

           **Arguments:**

           factor
                A scalar factor.
        '''
        check_type('factor', factor, float, int)
        self._array *= factor

    def get_element(self, i, j, k):
        '''Return a matrix element'''
        return self._array[i, j, k]

    def set_element(self, i, j, k, value):
        '''Set a matrix element'''
        self._array[i, j, k] = value

    #
    # Properties
    #

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._array.shape[0]

    nbasis = property(_get_nbasis)

    def _get_shape(self):
        '''The shape of the object'''
        return self._array.shape

    shape = property(_get_shape)

    #
    # New methods for this implementation
    # TODO: consider adding these to base class
    #

    def contract_two_to_two(self, subscripts, inp, out=None, factor=1.0, clear=True):
        '''Contracts self with two-index to obtain two-index.

           **Arguments:**

           subscripts
                One of ``abc,ab->ac``, ``abc,ab->ca``, ``abc,bc->ba``,
                ``abc,bc->ab``, ``abc,cb->ac``.

           inp
                A DenseTwoIndex input object.

           **Optional arguments:**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`
        '''
        check_options('subscripts', subscripts, 'abc,ab->ac', 'abc,ab->ca',
            'abc,bc->ba', 'abc,bc->ab', 'abc,cb->ac')
        return DenseLinalgFactory.einsum(subscripts, out, factor, clear, self, inp)

    # FIXME: can this be avoided by reorganizing with later contractions? This eats memory.
    def iadd_expand_two_one(self, subscripts, two, one, factor=1.0):
        '''In-place addition of expanded two-index and one-index to three-index.

           **Arguments:**

           subscripts
                Contraction type: ``ab,c->cab`` or ``ab,c->acb``.

           two
                A DenseTwoIndex object.

           one
                A DenseOneIndex object.

           **Optional arguments:**

           factor
                The added term is scaled by this factor.
        '''
        check_options('subscripts', subscripts, 'ab,c->cab', 'ab,c->acb')
        DenseLinalgFactory.einsum(subscripts, self, factor, False, two, one)

    def iadd_contract_three_two(self, subscripts, three, two, factor=1.0):
        '''In-place addition of a contraction of three-index with two-index.

           **Arguments:**

           subscripts
                Contraction type: ``abc,db->adc``, ``abc,dc->adb``,
                ``abc,db->dac``, ``abc,dc->dab``.

           three
                A DenseThreeIndex object.

           two
                A DenseTwoIndex object.

           **Optional arguments:**

           factor
                The added term is scaled by this factor.
        '''
        check_options('subscripts', subscripts, 'abc,db->adc', 'abc,dc->adb',
            'abc,db->dac', 'abc,dc->dab')
        DenseLinalgFactory.einsum(subscripts, self, factor, False, three, two)


class DenseFourIndex(FourIndex):
    """Dense symmetric four-dimensional matrix.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """

    #
    # Constructor and destructor
    #

    def __init__(self, nbasis, symmetry=8):
        """
           **Arguments:**

           nbasis
                The number of basis functions.

           **Optional arguments:**

           symmetry
                The supposed symmetry number of the four-index object.

                * ``8``: all possible symmtries with respect to index
                         permutations apply.
                * ``4``: symmetric with respect to swapping two indexes of the
                         same electron.
                * ``2``: symmetric with respect to swapping two electrons.
                * ``1``: no symmetries.
        """
        check_options('symmetry', symmetry, 1, 2, 4, 8)
        self._array = np.zeros((nbasis, nbasis, nbasis, nbasis))
        self._symmetry = symmetry
        log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)

    #
    # Methods from base class
    #

    def __check_init_args__(self, nbasis, symmetry=8):
        '''Is self compatible with the given constructor arguments?'''
        assert nbasis == self.nbasis
        assert symmetry == self.symmetry

    def __eq__(self, other):
        '''Compare self with other'''
        return isinstance(other, DenseFourIndex) and \
            other.nbasis == self.nbasis and \
            other.symmetry == self.symmetry and \
            (other._array == self._array).all()

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct an instance from data previously stored in an h5py.Group.

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        nbasis = grp['array'].shape[0]
        symmetry = grp.attrs['symmetry']
        result = cls(nbasis, symmetry)
        grp['array'].read_direct(result._array)
        return result

    def to_hdf5(self, grp):
        '''Dump this object in an h5py.Group

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        grp.attrs['class'] = self.__class__.__name__
        grp.attrs['symmetry'] = self.symmetry
        grp['array'] = self._array

    def new(self):
        '''Return a new four-index object with the same nbasis'''
        return DenseFourIndex(self.nbasis, self.symmetry)

    def _check_new_init_args(self, other):
        '''Check whether an already initialized object is compatible'''
        other.__check_init_args__(self.nbasis, self.symmetry)

    new.__check_init_args__ = _check_new_init_args

    def clear(self):
        '''Reset all elements to zero.'''
        self._array[:] = 0.0

    def copy(self, begin=0, end=None):
        '''Return a copy of (a part of) the object

           **Optional arguments:**

           begin, end
                Can be used to select a subblock of the object. When not given,
                the full range is used. The same begin and end is applied to
                all indexes to maintain the eight-fold symmetry.
        '''
        if end is None:
            end = self._array.shape[0]
        result = DenseFourIndex(end - begin)
        result._array[:] = self._array[begin:end, begin:end, begin:end, begin:end]
        return result

    def assign(self, other):
        '''Assign with the contents of another object

           **Arguments:**

           other
                Another DenseFourIndex object or a scalar.
        '''
        check_type('other', other, DenseFourIndex)
        self._array[:] = other._array
        self._symmetry = other.symmetry

    def randomize(self):
        '''Fill with random normal data'''
        self._array[:] = np.random.normal(0, 1, self.shape)
        if self.symmetry in (2, 8):
            self._array[:] += self._array.transpose(1,0,3,2)
            self._array *= 0.5
        if self.symmetry in (4, 8):
            self._array[:] += self._array.transpose(2,3,0,1)
            self._array *= 0.5
            self._array[:] += self._array.transpose(1,2,3,0)
            self._array *= 0.5

    def permute_basis(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.

           **Arguments:**

           permutation
                An integer numpy array that defines the new order of the basis
                functions.
        '''
        self._array[:] = self._array[permutation]
        self._array[:] = self._array[:,permutation]
        self._array[:] = self._array[:,:,permutation]
        self._array[:] = self._array[:,:,:,permutation]

    def change_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.

           **Arguments:**

           signs
                A numpy array with sign changes indicated by +1 and -1.
        '''
        self._array *= signs
        self._array *= signs.reshape(-1,1)
        self._array *= signs.reshape(-1,1,1)
        self._array *= signs.reshape(-1,1,1,1)

    def iadd(self, other, factor=1.0):
        '''Add another DenseFourIndex object in-place, multiplied by factor

           **Arguments:**

           other
                A DenseFourIndex instance to be added

           **Optional arguments:**

           factor
                The added term is scaled by this factor.
        '''
        check_type('other', other, DenseFourIndex)
        check_type('factor', factor, float, int)
        def from_mask(mask):
            return {0: 1, 1: 2, 2: 4, 3: 8}[mask]
        def to_mask(sym):
            return {1: 0, 2: 1, 4: 2, 8: 3}[sym]
        self._symmetry = from_mask(to_mask(self.symmetry) & to_mask(other.symmetry))
        self._array += other._array*factor

    def iscale(self, factor):
        '''In-place multiplication with a scalar

           **Arguments:**

           factor
                A scalar factor.
        '''
        check_type('factor', factor, float, int)
        self._array *= factor

    def get_element(self, i, j, k, l):
        '''Return a matrix element'''
        return self._array[i, j, k, l]

    def set_element(self, i, j, k, l, value):
        '''Set a matrix element

           If self.symmetry is larger than 1, other matrix elements are also
           set to maintain the proper symmetry.
        '''
        self._array[i,j,k,l] = value
        if self.symmetry in (2, 8):
            self._array[j,i,l,k] = value
        if self.symmetry in (4, 8):
            self._array[k,j,i,l] = value
            self._array[i,l,k,j] = value
            self._array[k,l,i,j] = value
        if self.symmetry == 8:
            self._array[l,k,j,i] = value
            self._array[j,k,l,i] = value
            self._array[l,i,j,k] = value

    #
    # Properties
    #

    def _get_nbasis(self):
        '''The number of basis functions'''
        return self._array.shape[0]

    nbasis = property(_get_nbasis)

    def _get_shape(self):
        '''The shape of the object'''
        return self._array.shape

    shape = property(_get_shape)

    def _get_symmetry(self):
        return self._symmetry

    symmetry = property(_get_symmetry)

    #
    # New methods for this implementation
    # TODO: consider adding these to base class
    #

    def check_symmetry(self):
        """Check the symmetry of the array."""
        if self.symmetry in (2, 8):
            assert (self._array == self._array.transpose(1,0,3,2)).all()
        if self.symmetry in (4, 8):
            assert (self._array == self._array.transpose(2,3,0,1)).all()
            assert (self._array == self._array.transpose(2,1,0,3)).all()
            assert (self._array == self._array.transpose(0,3,2,1)).all()
        if self.symmetry == 8:
            assert (self._array == self._array.transpose(3,2,1,0)).all()
            assert (self._array == self._array.transpose(3,0,1,2)).all()
            assert (self._array == self._array.transpose(1,2,3,0)).all()

    def itranspose(self):
        '''In-place transpose: ``0,1,2,3 -> 1,0,3,2``'''
        self._array[:] = self._array.transpose(1,0,3,2)

    def sum(self):
        '''Return the sum of all elements'''
        return np.sum(self._array)

    def iadd_exchange(self):
        '''In-place addition of its own exchange contribution'''
        self._array -= np.einsum('abcd->abdc', self._array)
        self._symmetry = {8:2, 4:1, 2:2, 1:1}[self._symmetry]

    def slice_to_two(self, subscripts, out=None, factor=1.0, clear=True):
        """Returns a two-index contraction of the four-index object.

           **Arguments:**

           superscripts
                Any of ``aabb->ab``, ``abab->ab``, ``abba->ab``

           **Optional arguments:**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`
        """
        check_options('subscripts', subscripts, 'aabb->ab', 'abab->ab', 'abba->ab')
        return DenseLinalgFactory.einsum(subscripts, out, factor, clear, self)

    def slice_to_three(self, subscripts, out=None, factor=1.0, clear=True):
        """Returns a three-index contraction of the four-index object.

           **Arguments:**

           superscripts
                Any of ``abcc->bac``, ``abcc->abc``, ``abcb->abc``, ``abbc->abc``

           **Optional arguments:**

           out, factor, clear
                See :py:meth:`DenseLinalgFactory.einsum`
        """
        check_options('subscripts', subscripts, 'abcc->bac', 'abcc->abc', 'abcb->abc', 'abbc->abc')
        return DenseLinalgFactory.einsum(subscripts, out, factor, clear, self)

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
        return DenseLinalgFactory.einsum(subscripts, out, factor, clear, self, two)

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
#        return DenseLinalgFactory.einsum(subscripts, out, factor, clear, self, two)
        if subscripts == 'abcd,bd->ac':
            return DenseLinalgFactory.tensordot(self, two, ([1,3], [1,0]), out, factor, clear)
        elif subscripts == 'abcd,cb->ad':
            return DenseLinalgFactory.tensordot(self, two, ([1,2], [1,0]), out, factor, clear)

    def assign_four_index_transform(self, ao_integrals, exp0, exp1=None, exp2=None, exp3=None, method='tensordot'):
        '''Perform four index transformation.

           **Arguments:**

           oa_integrals
                A DenseFourIndex with integrals in atomic orbitals.

           exp0
                A DenseExpansion object with molecular orbitals

           **Optional arguments:**

           exp1, exp2, exp3
                Can be provided to transform each index differently.

           method
                Either ``einsum`` or ``tensordot`` (default).
        '''
        # parse arguments
        check_type('ao_integrals', ao_integrals, DenseFourIndex)
        exp0, exp1, exp2, exp3 = parse_four_index_transform_exps(exp0, exp1, exp2, exp3, DenseExpansion)
        # update self.symmetry
        if not (exp0 is exp1 and exp2 is exp3):
            # 8->4, 2->1
            self._symmetry = {8:4, 4:4, 2:1, 1:1}[self.symmetry]
        if not (exp0 is exp2 and exp1 is exp3):
            # 8->2, 4->1
            self._symmetry = {8:2, 4:1, 2:2, 1:1}[self.symmetry]
        # actual transform
        if method == 'einsum':
            # The order of the dot products is according to literature
            # conventions.
            self._array[:] = np.einsum('sd,pqrs->pqrd', exp3.coeffs, ao_integrals._array, casting='no', order='C')
            self._array[:] = np.einsum('rc,pqrd->pqcd', exp2.coeffs, self._array, casting='no', order='C')
            self._array[:] = np.einsum('qb,pqcd->pbcd', exp1.coeffs, self._array, casting='no', order='C')
            self._array[:] = np.einsum('pa,pbcd->abcd', exp0.coeffs, self._array, casting='no', order='C')
        elif method == 'tensordot':
            # because the way tensordot works, the order of the dot products is
            # not according to literature conventions.
            self._array[:] = np.tensordot(ao_integrals._array, exp0.coeffs, axes=([0],[0]))
            self._array[:] = np.tensordot(self._array, exp1.coeffs, axes=([0],[0]))
            self._array[:] = np.tensordot(self._array, exp2.coeffs, axes=([0],[0]))
            self._array[:] = np.tensordot(self._array, exp3.coeffs, axes=([0],[0]))
        else:
            raise ValueError('The method must either be \'einsum\' or \'tensordot\'.')
