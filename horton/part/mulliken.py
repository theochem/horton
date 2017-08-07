# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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
"""Mulliken partitioning"""


import numpy as np


__all__ = ["partition_mulliken", "get_mulliken_operators"]


def get_shell_nbasis(shell_type):
    """Return number of basis functions in the given shell type.

    Parameters
    ----------
    shell_type : int
        Integer representing the shell.
    """
    if shell_type > 0:
        # cartesian
        return (shell_type + 1) * (shell_type + 2) / 2
    elif shell_type == -1:
        # should not happen
        return -1
    else:
        # pure
        return -2 * shell_type + 1


def partition_mulliken(operator, nbasis, shell_types, shell_maps, index):
    """Fill in the mulliken operator in the first argument.

    Parameters
    ----------
    operator : np.ndarray, shape=(nbasis, nbasis), dtype=float
        A Two index operator to which the Mulliken mask is applied.
    nbasis : int
        Number of basis function.
    shell_types : list
        Sequence of integers representing the basis shell types.
    shell_maps : list
        Sequence of integers representing the center each shell belongs to.
    index : int
        The index of the atom (center) for which the Mulliken operator needs to be
        constructed.

    This routine implies that the first ``natom`` centers in the obasis corresponds to the
    atoms in the system.
    """
    mask = np.zeros(nbasis, dtype=bool)
    begin = 0
    for ishell, shell_type in enumerate(shell_types):
        end = begin + get_shell_nbasis(shell_type)
        # check whether shell belongs to center denoted by index
        if shell_maps[ishell] != index:
            mask[begin:end] = True
        begin = end
    operator[mask] = 0.0
    operator[:] = 0.5 * (operator + operator.T)


def get_mulliken_operators(overlap, ncenter, shell_types, shell_maps):
    """Return a list of Mulliken operators for the given basis.

    Parameters
    ----------
    overlap : np.ndarray, shape=(nbasis, nbasis), dtype=float
        The overlap matrix in a given basis.
    ncenter : int
        Number of basis centers.
    shell_types : list
        Sequence of integers representing the basis shell types.
    shell_maps : list
        Sequence of integers representing the center each shell belongs to.
    """
    nbasis = len(overlap)
    operators = []
    for icenter in range(ncenter):
        operator = overlap.copy()
        partition_mulliken(operator, nbasis, shell_types, shell_maps, icenter)
        operators.append(operator)
    return operators
