# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
"""Base classes"""


from horton.utils import check_type


__all__ = [
    'LinalgFactory', 'LinalgObject', 'OneIndex', 'Expansion', 'TwoIndex',
    'ThreeIndex', 'FourIndex',
    'parse_four_index_transform_exps',
]


# Possible symmetries for four-index objects

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

    @classmethod
    def from_hdf5(cls, grp):
        '''Construct an instance from data previously stored in an h5py.Group.

           **Arguments:**

           grp
                An h5py.Group object.
        '''
        default_nbasis = grp.attrs.get('default_nbasis')
        return cls(default_nbasis)

    def to_hdf5(self, grp):
        '''Write a LinalgFactory to an HDF5 group

           **Argument:**

           grp
                A h5py.Group instance to write to.
        '''
        grp.attrs['class'] = self.__class__.__name__
        if self.default_nbasis is not None:
            grp.attrs['default_nbasis'] = self.default_nbasis

    def set_default_nbasis(self, nbasis):
        self.default_nbasis = nbasis

    def create_one_index(self, nbasis=None):
        raise NotImplementedError

    def create_expansion(self, nbasis=None):
        raise NotImplementedError

    def create_two_index(self, nbasis=None):
        raise NotImplementedError

    def create_four_index(self, nbasis=None):
        raise NotImplementedError

    def create_three_index(self, nbasis=None):
        raise NotImplementedError


class LinalgObject(object):
    def __check_init_args__(self, nbasis):
        raise NotImplementedError

    def __eq__(self, other):
        raise NotImplementedError

    @classmethod
    def from_hdf5(cls, grp):
        raise NotImplementedError

    def to_hdf5(self, grp):
        raise NotImplementedError

    def new(self):
        raise NotImplementedError

    def __clear__(self):
        '''Part of the API specified in horton.cache'''
        self.clear()

    def clear(self):
        raise NotImplementedError

    def copy(self):
        raise NotImplementedError

    def assign(self, other):
        raise NotImplementedError

    def randomize(self):
        raise NotImplementedError

    def permute_basis(self, permutation):
        raise NotImplementedError

    def change_basis_signs(self, signs):
        raise NotImplementedError


class NIndexObject(LinalgObject):
    def __init__(self, nbasis):
        raise NotImplementedError

    def iadd(self, other, factor=1):
        raise NotImplementedError

    def iscale(self, factor):
        raise NotImplementedError

    def _get_shape(self):
        raise NotImplementedError

    shape = property(_get_shape)

    def _get_ndim(self):
        '''The number of axes in the N-index object.'''
        return len(self.shape)

    ndim = property(_get_ndim)

    def _fix_ends(self, *ends):
        shape = self.shape
        if len(shape) != len(ends):
            raise TypeError('The argument \'ends\' must have the same length as \'self.shape\'.')
        return tuple([ends[i] or shape[i] for i in xrange(len(shape))])


class OneIndex(NIndexObject):
    def copy(self, begin=0, end=None):
        raise NotImplementedError

    def get_element(self, i):
        raise NotImplementedError

    def set_element(self, i, value):
        raise NotImplementedError


class Expansion(LinalgObject):
    def __init__(self, nbasis, nfn=None):
        raise NotImplementedError

    def check_normalization(self, olp, eps=1e-4):
        raise NotImplementedError

    def error_eigen(self, fock, overlap):
        raise NotImplementedError


class TwoIndex(NIndexObject):
    def __init__(self, nbasis):
        raise NotImplementedError

    def copy(self, begin0=0, end0=None, begin1=None, end1=None):
        raise NotImplementedError

    def get_element(self, i, j):
        raise NotImplementedError

    def set_element(self, i, j, value):
        raise NotImplementedError

    def sum(self):
        raise NotImplementedError

    def trace(self):
        raise NotImplementedError

    def itranspose(self):
        raise NotImplementedError

    def inner(self, vec0, vec1):
        raise NotImplementedError


class ThreeIndex(NIndexObject):
    def copy(self, begin0=0, end0=None, begin1=None, end1=None, begin2=None, end2=None):
        raise NotImplementedError

    def get_element(self, i, j, k):
        raise NotImplementedError

    def set_element(self, i, j, k, value):
        raise NotImplementedError


class FourIndex(NIndexObject):
    def copy(self, begin=0, end=None):
        raise NotImplementedError

    def get_element(self, i, j, k, l):
        raise NotImplementedError

    def set_element(self, i, j, k, l, value):
        raise NotImplementedError


def parse_four_index_transform_exps(exp0, exp1, exp2, exp3, Class):
    '''Parse the optional arguments exp1, exp2 and exp3.

       **Arguments:**

       exp0, exp1, exp2, exp3
            Four sets of orbitals for the mo transformation. Some may be None
            but only the following not None combinations are allowed:

            * ``(exp0,)``: maintain eight-fold symmetry (if any)
            * ``(exp0, exp1)``: maintain four-fold symmetry (if any)
            * ``(exp0, exp2)``: maintain two-fold symmetry (if any)
            * ``(exp0, exp1, exp2, exp3)``: break all symmetry

       Class
            The expected class of the exps objects.


       **Returns:** exp0, exp1, exp2, exp3. (All not None)
    '''
    # Four supported situations
    if exp1 is None and exp2 is None and exp3 is None:
        # maintains eight-fold symmetry
        exp1 = exp0
        exp2 = exp0
        exp3 = exp0
    elif exp2 is None and exp3 is None:
        # maintains four-fold symmetry
        exp2 = exp0
        exp3 = exp1
    elif exp1 is None and exp3 is None:
        # maintains two-fold symmetry
        exp1 = exp0
        exp3 = exp2
    elif exp1 is None or exp2 is None or exp3 is None:
        # the only other allowed case is no symmetry.
        raise TypeError('It is not clear how to interpret the optional arguments exp1, exp2 and exp3.')
    check_type('exp0', exp0, Class)
    check_type('exp1', exp1, Class)
    check_type('exp2', exp2, Class)
    check_type('exp3', exp3, Class)
    return exp0, exp1, exp2, exp3
