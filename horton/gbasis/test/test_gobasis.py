# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from horton import *


def test_shell_nbasis():
    assert get_shell_nbasis(-3) == 7
    assert get_shell_nbasis(-2) == 5
    assert get_shell_nbasis( 0) == 1
    assert get_shell_nbasis( 1) == 3
    assert get_shell_nbasis( 2) == 6
    assert get_shell_nbasis( 3) == 10
    try:
        get_shell_nbasis(-1)
        assert False
    except ValueError:
        pass


def test_gobasis_consistency():
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 0, 1, 1, 1, 1])
    nprims = np.array([2, 3, 3, 5, 5, 5, 7])
    shell_types = np.array([2, 1, 0, -2, 3, 0, 1])
    alphas = np.random.uniform(0, 1, nprims.sum())
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())

    gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    assert gobasis.nbasis == 29
    assert gobasis.max_shell_type == 3
    scales = gobasis.get_scales()
    assert abs(scales[0] - gob_normalization(alphas[0], np.array([2, 0, 0]))) < 1e-10

    shell_types = np.array([1, 1, 0, -2, -2, 0, 1])
    gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    assert gobasis.nbasis == 21
    assert gobasis.max_shell_type == 2

    # The center indexes in the shell_map are out of range.
    shell_map[0] = 2
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    shell_map[0] = 0

    # The size of the array shell_types does not match the sum of nprims.
    shell_types = np.array([1, 1])
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except TypeError:
        pass
    shell_types = np.array([1, 1, 0, -2, -2, 0, 1])

    # The elements of nprims should be at least 1.
    nprims[1] = 0
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    nprims[1] = 3

    # The size of the array alphas does not match the sum of nprims.
    alphas = np.random.uniform(-1, 1, 2)
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except TypeError:
        pass
    alphas = np.random.uniform(-1, 1, nprims.sum())

    # Encountered the nonexistent shell_type -1.
    shell_types[1] = -1
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    shell_types[1] = 1

    # The size of con_coeffs does not match nprims.
    con_coeffs = np.random.uniform(-1, 1, 3)
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except TypeError:
        pass
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())

    # Exceeding the maximym shell type (above):
    shell_types[0] = get_max_shell_type()+1
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    shell_types[0] = 2

    # Exceeding the maximym shell type (below):
    shell_types[0] = -get_max_shell_type()-1
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    shell_types[0] = 2
