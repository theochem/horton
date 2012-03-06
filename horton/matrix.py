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
"""


import numpy as np


class Dense2(object):
    """Dense symmetric two-dimensional matrix.

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


class Dense4(object):
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
