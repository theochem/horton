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
"""Cholesky decomposition of four-index objects"""


import numpy as np
from horton.log import log
from horton.matrix.cext import compute_slice_abcc, compute_slice_abbc, subtract_slice_abbc
from horton.matrix.dense import DenseLinalgFactory, DenseTwoIndex, DenseFourIndex


__all__ = [
    'CholeskyFourIndex', 'CholeskyLinalgFactory',
]


class CholeskyLinalgFactory(DenseLinalgFactory):
    def create_four_index(self, nbasis=None, nvec=None, array=None):
        nbasis = nbasis or self.default_nbasis

        if array is not None:
            if nvec is not None:
                assert array.shape[0] == nvec
            return CholeskyFourIndex(nbasis=nbasis, array=array)
        elif nvec is not None:
            return CholeskyFourIndex(nbasis=nbasis, nvec=nvec)
        else:
            raise NotImplementedError

    def _check_four_index_init_args(self, four_index, nbasis=None, nvec=None, array=None):
        nbasis = nbasis or self.default_nbasis
        four_index.__check_init_args__(nbasis)

    create_four_index.__check_init_args__ = _check_four_index_init_args


class CholeskyFourIndex(DenseFourIndex):
    """Dense symmetric four-dimensional matrix.

       This is the most inefficient implementation in terms of memory usage and
       computer time. Due to its simplicity, it is trivial to implement. This
       implementation mainly serves as a reference for testing purposes.
    """
    def __init__(self, nbasis, nvec=None, array=None):
        """
           **Arguments:**

           nbasis
                The number of basis functions.
        """
        self._array = None
        self._array2 = self._array

        if array is not None:
            self.assign_array(array)
        elif nvec is not None:
            self.assign_array(np.zeros([nvec, nbasis, nbasis]))

        if self._array is not None:
            log.mem.announce(self._array.nbytes)

    def __del__(self):
        if log is not None:
            if hasattr(self, '_array'):
                log.mem.denounce(self._array.nbytes)
                if self._array2 is not self._array:
                    log.mem.denounce(self._array2.nbytes)

    def __check_init_args__(self, nbasis):
        assert nbasis == self.nbasis
        assert self._array is not None

    def reset_array2(self):
        """ Deallocates the second cholesky vector and sets it to match the first.
        """
        self._array2 = self._array

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
        return self._array.shape[1]

    nbasis = property(_get_nbasis)

    def set_element(self, i, j, k, l, value):
        #    <ij|kl> = <ji|lk> = <kl|ij> = <lk|ji> =
        #    <il|kj> = <jk|li> = <kj|il> = <li|jk>
        raise NotImplementedError

    def get_element(self, i, j, k, l):
        return np.dot(self._array[:,i,k], self._array2[:,j,l])

    def get_slice(self, indices, out=None, subtract=False, clear=True):
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

            ** Optional Arguments **
                out
                    A numpy array for storing the results.
                factor
                    A scaling factor to multiply the results being stored to out
                clear
                    A boolean that determines whether the output array is
                    cleared. Defaults to True.
        """

        indc_in,indc_out = "".join(indices.split()).split("->")
        assert len(indc_in) == 4
        assert set(indc_in) == set(indc_out)
        assert len(set("xyz") & set(indc_in)) == 0


        if out is None:
            dims = [self.nbasis, self.nbasis]
            if len(indc_out) == 3:
                dims += [self.nbasis]
            out = np.zeros(dims)
        elif clear:
            out[:] = 0.0

        if indices == 'abbc->abc': #abbc
            if subtract:
                subtract_slice_abbc(self._array, self._array2, out, self.nbasis,
                                        self._array.shape[0])
            else:
                compute_slice_abbc(self._array, self._array2, out, self.nbasis,
                                        self._array.shape[0])
#            for i in np.arange(self._array.shape[0]):
#                out += np.einsum('ab,bc->abc',self._array[i,:,:],
#                                        self._array2[i,:,:])
        elif subtract:
            raise NotImplementedError
        elif indices == 'abcc->bac': #acbc
            print self._array.shape, self._array2.shape, out.shape, self.nbasis
            compute_slice_abcc(self._array, self._array2, out, self.nbasis,
                                                    self._array.shape[0])
            out = out.swapaxes(0,1) #DOES THIS WORK?
#            for i in np.arange(self._array.shape[0]):
#                out += np.einsum('ac,bc->bac', self._array[i,:,:],
#                                        self._array2[i,:,:])
        elif indices == 'abcc->abc': #acbc
            compute_slice_abcc(self._array, self._array2, out, self.nbasis,
                                self._array.shape[0])
#            for i in np.arange(self._array.shape[0]):
#                out += np.einsum('ac,bc->abc', self._array[i,:,:],
#                                        self._array2[i,:,:])
        elif indices == 'abcb->abc': #acbb
            L_r = np.diagonal(self._array2, axis1=1, axis2=2)
            out[:] = np.tensordot(self._array, L_r, [0,0]).swapaxes(1,2)
        else:
            #error checking
            idx_string = "x{},x{}->{}".format(indc_in[0]+indc_in[2],
                                               indc_in[1]+indc_in[3],
                                               indc_out)

            out[:] = np.einsum(idx_string, self._array, self._array2)
        return out

    def _get_dense(self):
        ''' ONLY FOR TESTING. Super expensive operation!
        '''
        return np.einsum('kac,kbd->abcd', self._array, self._array2)

    def assign(self, other):
        if not isinstance(other, CholeskyFourIndex):
            raise TypeError('The other object must also be CholeskyFourIndex instance. . Got ', type(other), ' instead.')
        if self._array is None:
            self._array = np.ndarray(other._array.shape)
            log.mem.announce(self._array.nbytes)

        self._array[:] = other._array
        if other._array is other._array2:
            #arrays are the same
            self._array2 = self._array
        else:
            #arrays have been transformed
            if self._array is None:
                self._array2 = np.ndarray(other._array2.shape)
                log.mem.announce(self._array2.nbytes)

            self._array2[:] = other._array2

    def assign_array(self, other):
        ''''''
        if not isinstance(other, np.ndarray):
            raise TypeError('The other object must be np.ndarray instance. . Got ', type(other), ' instead.')

        self._array = other.copy() #maybe copy not needed?
        self._array2 = self._array

    def copy(self):
        '''Return a copy of the current four-index operator'''
        result = CholeskyFourIndex(self.nbasis)
        result.assign(self)
        return result

    def iscale(self, factor):
        self._array *= np.sqrt(factor)

        if self._array is not self._array2:
            #arrays have been transformed
            self._array2 *= np.sqrt(factor)

    def check_symmetry(self):
        """Check the symmetry of the array."""
        raise NotImplementedError

    def transpose(self):
        if self._array is not self._array2:
            #arrays have been transformed
            temp = self._array
            self._array = self._array2
            self._array2 = temp

    def esum(self):
        return np.tensordot(self._array, self._array2,(0,0)).sum() #expensive!!

    def apply_direct(self, dm, output):
        """Compute the direct dot product with a density matrix."""
        if not isinstance(dm, DenseTwoIndex):
            raise TypeError('The dm argument must be a DenseTwoIndex class. Got ', type(dm), ' instead.')
        if not isinstance(output, DenseTwoIndex):
            raise TypeError('The output argument must be a DenseTwoIndex class. Got ', type(output), ' instead.')
        result = np.tensordot(self._array, np.tensordot(self._array2, dm._get_dense(), axes=([(1,2),(1,0)])), [0,0])
        output.assign_array(result)

    def apply_exchange(self, dm, output):
        """Compute the exchange dot product with a density matrix."""
        if not isinstance(dm, DenseTwoIndex):
            raise TypeError('The dm argument must be a DenseTwoIndex class. Got ', type(dm), ' instead.')
        if not isinstance(output, DenseTwoIndex):
            raise TypeError('The output argument must be a DenseTwoIndex class. Got ', type(output), ' instead.')
        result = np.tensordot(self._array, np.tensordot(self._array2, dm._get_dense(), axes=([2,1])), ([0,2],[0,2]))
        output.assign_array(result)

    def apply_four_index_transform_tensordot(self, ao_integrals, aorb, aorb2=None, aorb3=None, aorb4=None):
        '''Perform four index transformation using np.tensordot.
        '''
        if aorb2 is None and aorb4 is None:
            aorb2 = aorb
        else:
            aorb2 = aorb4

        if aorb3 is None:
            aorb3 = aorb
        if aorb4 is None:
            aorb4 = aorb2

        if aorb != aorb3 or aorb2 != aorb4:
            raise NotImplementedError

        if not isinstance(ao_integrals, CholeskyFourIndex):
            raise TypeError('The AO integral argument must be a CholeskyFourIndex class. Got ', type(ao_integrals), ' instead.')
        self._array[:] = np.tensordot(ao_integrals._array, aorb2.coeffs, axes=([1],[0]))
        self._array[:] = np.tensordot(self._array, aorb2.coeffs, axes=([1],[0]))

        if aorb != aorb2 and aorb3 != aorb4:
            if self._array is self._array2:
                #must allocate memory first
                self._array2 = np.zeros_like(self._array)
                log.mem.announce(self._array2.nbytes)
            self._array2[:] = np.tensordot(ao_integrals._array, aorb.coeffs, axes=([1],[0]))
            self._array2[:] = np.tensordot(self._array2, aorb.coeffs, axes=([1],[0]))


    def apply_four_index_transform_einsum(self, ao_integrals, aorb, aorb2=None, aorb3=None, aorb4=None):
        '''Perform four index transformation using np.einsum.
        '''
        if aorb2 is None and aorb4 is None:
            aorb2 = aorb
        else:
            aorb2 = aorb4

        if aorb3 is None:
            aorb3 = aorb
        if aorb4 is None:
            aorb4 = aorb2

        if aorb != aorb3 or aorb2 != aorb4:
            raise NotImplementedError

        if not isinstance(ao_integrals, CholeskyFourIndex):
            raise TypeError('The AO integral argument must be a CholeskyFourIndex class. Got ', type(ao_integrals), ' instead.')
        self._array[:] = np.einsum('ai,kab->kib',aorb2.coeffs,ao_integrals._array)
        self._array[:] = np.einsum('bj,kib->kij',aorb2.coeffs,self._array)

        if aorb != aorb2 and aorb3 != aorb4:
            if self._array is self._array2:
                #must allocate memory first
                self._array2 = np.zeros_like(self._array)
                log.mem.announce(self._array2.nbytes)
            self._array2[:] = np.einsum('ai,kab->kib',aorb.coeffs,ao_integrals._array2)
            self._array2[:] = np.einsum('bj,kib->kij',aorb.coeffs,self._array2)

    def clear(self):
        self._array[:] = 0.0
        if self._array is not self._array2:
            self._array2[:] = 0.0

    def apply_basis_permutation(self, permutation):
        '''Reorder the coefficients for a given permutation of basis functions.
        '''
        raise NotImplementedError

    def apply_basis_signs(self, signs):
        '''Correct for different sign conventions of the basis functions.'''
        raise NotImplementedError

    def add_exchange_part(self):
        '''Sort four-index exchange integrals for OAP1roG (<ij||kj> -> <ijk>)
        '''
        raise NotImplementedError
