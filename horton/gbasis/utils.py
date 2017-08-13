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
"""Utility Functions"""


from os import path
import numpy as np


def typecheck_geo(coordinates=None, numbers=None, pseudo_numbers=None,
                  need_coordinates=True, need_numbers=True,
                  need_pseudo_numbers=True):
    """Type check a molecular geometry specification

       **Arguments:**

       coordinates
            A (N, 3) float array with Cartesian coordinates of the atoms.

       numbers
            A (N,) int vector with the atomic numbers.

       **Optional arguments:**

       pseudo_numbers
            A (N,) float array with pseudo-potential core charges.

       need_coordinates
            When set to False, the coordinates can be None, are not type checked
            and not returned.

       need_numbers
            When set to False, the numbers can be None, are not type checked
            and not returned.

       need_pseudo_numbers
            When set to False, the pseudo_numbers can be None, are not type
            checked and not returned.

       **Returns:** ``[natom]`` + all arguments that were type checked. The
       pseudo_numbers argument is converted to a floating point array.
    """
    # Determine natom
    if coordinates is not None:
        natom = len(coordinates)
    elif numbers is not None:
        natom = len(numbers)
    elif pseudo_numbers is not None:
        natom = len(pseudo_numbers)
    else:
        raise TypeError('At least one argument is required and should not be None')

    # Typecheck coordinates:
    if coordinates is None:
        if need_coordinates:
            raise TypeError('Coordinates can not be None.')
    else:
        if coordinates.shape != (natom, 3) or not np.issubdtype(coordinates.dtype, np.float):
            raise TypeError('The argument centers must be a float array with shape (natom,3).')

    # Typecheck numbers
    if numbers is None:
        if need_numbers:
            raise TypeError('Numbers can not be None.')
    else:
        if numbers.shape != (natom,) or not np.issubdtype(numbers.dtype, np.integer):
            raise TypeError('The argument numbers must be a vector with length natom.')

    # Typecheck pseudo_numbers
    if pseudo_numbers is None:
        if need_pseudo_numbers:
            pseudo_numbers = numbers.astype(float)
    else:
        if pseudo_numbers.shape != (natom,):
            raise TypeError('The argument pseudo_numbers must be a vector with length natom.')
        if not np.issubdtype(pseudo_numbers.dtype, np.float):
            pseudo_numbers = pseudo_numbers.astype(float)

    # Collect return values
    result = [natom, ]
    if need_coordinates:
        result.append(coordinates)
    if need_numbers:
        result.append(numbers)
    if need_pseudo_numbers:
        result.append(pseudo_numbers)
    return result


def to_basis_path(fn):
    cur_pth = path.split(__file__)[0]
    pth = cur_pth + "/bsets/{}".format(fn)
    return pth
