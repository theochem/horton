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
"""Base classes"""


__all__ = [
    'LinalgFactory', 'LinalgObject', 'Expansion', 'OneIndex', 'TwoIndex',
    'ThreeIndex', 'FourIndex',
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

    def set_default_nbasis(self, nbasis):
        self.default_nbasis = nbasis

    def create_expansion(self, nbasis=None):
        raise NotImplementedError

    def create_two_index(self, nbasis=None, nfn=None):
        raise NotImplementedError

    def create_four_index(self, nbasis=None):
        raise NotImplementedError

    def create_three_index(self, nbasis=None):
        raise NotImplementedError

    def error_eigen(self, ham, overlap, expansion, epsilons):
        raise NotImplementedError

    def diagonalize(self, ham, overlap, expansion, epsilons):
        raise NotImplementedError

    def get_memory_two_index(self, nbasis=None):
        raise NotImplementedError

    def get_memory_four_index(self, nbasis=None):
        raise NotImplementedError


class LinalgObject(object):
    def __check_init_args__(self, nbasis):
        raise NotImplementedError

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
        '''Part of the API specified in horton.cache'''
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


class NIndexObject(LinalgObject):
    def __init__(self, nbasis):
        raise NotImplementedError

    def iadd(self, other, factor=1):
        raise NotImplementedError

    def iscale(self, factor):
        raise NotImplementedError


class OneIndex(NIndexObject):
    def set_element(self, i, value):
        raise NotImplementedError

    def get_element(self, i):
        raise NotImplementedError


class TwoIndex(NIndexObject):
    def set_element(self, i, j, value):
        raise NotImplementedError

    def get_element(self, i, j):
        raise NotImplementedError

    def expectation_value(self, dm):
        raise NotImplementedError

    def trace(self):
        raise NotImplementedError

    def itranspose(self):
        raise NotImplementedError

    def dot(self, vec0, vec1):
        raise NotImplementedError


class ThreeIndex(NIndexObject):
    def set_element(self, i, j, k, value):
        raise NotImplementedError

    def get_element(self, i, j, k):
        raise NotImplementedError


class FourIndex(NIndexObject):
    def set_element(self, i, j, k, l, value):
        raise NotImplementedError

    def get_element(self, i, j, k, l):
        raise NotImplementedError
