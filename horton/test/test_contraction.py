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



def get_test_i2gob():
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 1, 1])
    ncons = np.array([1, 2, 3, 1])
    nexps = np.array([2, 3, 5, 7])
    con_types = np.array([2, 1, 0, -2, 3, 0, 1])
    exponents = np.random.uniform(-1, 1, nexps.sum())
    con_coeffs = np.random.uniform(-1, 1, np.dot(nexps, ncons))
    nbasis = sum(get_con_nbasis(con_type) for con_type in con_types)
    return (
        I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis),
        centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis
    )


def test_i2gob_init():
    i2, centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis = get_test_i2gob()
    assert i2.max_nbasis == 10

    con_types = np.array([1, 1, 0, -2, -2, 0, 1])
    nbasis = sum(get_con_nbasis(con_type) for con_type in con_types)
    i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis)
    assert i2.max_nbasis == 6

    # -2) The center indexes in the shell_map are out of range.
    shell_map[0] = 2
    try:
        i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis)
        assert False
    except ValueError:
        pass
    shell_map[0] = 0

    # -3) The elements of ncons should be at least 1.
    ncons[3] = 0
    try:
        i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis)
        assert False
    except ValueError:
        pass
    ncons[3] = 3

    # -4) The size of the array con_types does not match the sum of ncons.
    con_types = np.array([1, 1])
    try:
        i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis)
        assert False
    except ValueError:
        pass
    con_types = np.array([1, 1, 0, -2, -2, 0, 1])

    # -5) The elements of nexps should be at least 1.
    nexps[1] = 0
    try:
        i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis)
        assert False
    except ValueError:
        pass
    nexps[1] = 3

    # -6) The size of the array exponents does not match the sum of nexps.
    exponents = np.random.uniform(-1, 1, 2)
    try:
        i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis)
        assert False
    except ValueError:
        pass
    exponents = np.random.uniform(-1, 1, nexps.sum())

    # -7) Encountered the nonexistent con_type -1.
    con_types[1] = -1
    try:
        i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis)
        assert False
    except ValueError:
        pass
    con_types[1] = 1

    # -8) The size of con_coeffs does not match nexp and ncon.
    con_coeffs = np.random.uniform(-1, 1, 3)
    try:
        i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, nbasis)
        assert False
    except ValueError:
        pass
    con_coeffs = np.random.uniform(-1, 1, np.dot(nexps, ncons))

    # -9) The total number of basis functions does not match the output array.
    try:
        i2 = I2Gob(centers, shell_map, ncons, nexps, con_types, exponents, con_coeffs, 5)
        assert False
    except ValueError:
        pass


def test_i2gob_inc_shell():
    i2 = get_test_i2gob()[0]
    assert i2.private_fields == (0, 0, 1, 1, 0, 0)
    assert i2.inc_shell() is True
    assert i2.private_fields == (0, 1, 1, 2, 0, 0)
    assert i2.inc_shell() is True
    assert i2.private_fields == (1, 1, 2, 2, 0, 0)
    assert i2.inc_shell() is True
    assert i2.private_fields == (0, 2, 1, 3, 0, 0)
    assert i2.inc_shell() is True
    assert i2.private_fields == (1, 2, 2, 3, 0, 0)
    assert i2.inc_shell() is True
    assert i2.private_fields == (2, 2, 3, 3, 0, 0)
    assert i2.inc_shell() is True
    assert i2.private_fields == (0, 3, 1, 1, 0, 0)


def test_i2gob_inc_con():
    i2 = get_test_i2gob()[0]
    assert i2.private_fields == (0, 0, 1, 1, 0, 0)
    assert i2.inc_con() is False
    assert i2.private_fields == (0, 0, 1, 1, 0, 0)
    assert i2.inc_shell() is True
    assert i2.private_fields == (0, 1, 1, 2, 0, 0)
    assert i2.inc_con() is True
    assert i2.private_fields == (0, 1, 1, 2, 0, 1)
    assert i2.inc_con() is False
    assert i2.private_fields == (0, 1, 1, 2, 0, 0)
    assert i2.inc_shell() is True
    assert i2.private_fields == (1, 1, 2, 2, 0, 0)
    assert i2.inc_con() is True
    assert i2.private_fields == (1, 1, 2, 2, 1, 0)
    assert i2.inc_con() is True
    assert i2.private_fields == (1, 1, 2, 2, 0, 1)
    assert i2.inc_con() is True
    assert i2.private_fields == (1, 1, 2, 2, 1, 1)


def test_i1pow_inc_l0():
    indexes = (0,0,0)
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,0,0)
    assert result is False


def test_i1pow_inc_l1():
    indexes = (1,0,0)
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,1,0)
    assert result is True
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,0,1)
    assert result is True
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (1,0,0)
    assert result is False


def test_i1pow_inc_l2():
    indexes = (2,0,0)
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (1,1,0)
    assert result is True
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (1,0,1)
    assert result is True
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,2,0)
    assert result is True
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,1,1)
    assert result is True
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (0,0,2)
    assert result is True
    indexes, result = i1pow_inc(*indexes)
    assert indexes == (2,0,0)
    assert result is False


def test_i2_pow_l0l0():
    i2p = I2Pow(0, 0, 3)
    assert i2p.fields == (0, 0, 0, 0, 0, 0, 0)
    assert i2p.inc() is False
    assert i2p.fields == (0, 0, 0, 0, 0, 0, 0)


def test_i2_pow_l1l0():
    i2p = I2Pow(1, 0, 3)
    assert i2p.fields == (1, 0, 0, 0, 0, 0, 0)
    assert i2p.inc() is True
    assert i2p.fields == (0, 1, 0, 0, 0, 0, 1)
    assert i2p.inc() is True
    assert i2p.fields == (0, 0, 1, 0, 0, 0, 2)
    assert i2p.inc() is False
    assert i2p.fields == (1, 0, 0, 0, 0, 0, 0)


def test_i2_pow_l2l1():
    i2p = I2Pow(2, 1, 10)
    assert i2p.fields == (2, 0, 0, 1, 0, 0, 0)
    assert i2p.inc() is True
    assert i2p.fields == (1, 1, 0, 1, 0, 0, 1)
    assert i2p.inc() is True
    assert i2p.fields == (1, 0, 1, 1, 0, 0, 2)
    assert i2p.inc() is True
    assert i2p.fields == (0, 2, 0, 1, 0, 0, 3)
    assert i2p.inc() is True
    assert i2p.fields == (0, 1, 1, 1, 0, 0, 4)
    assert i2p.inc() is True
    assert i2p.fields == (0, 0, 2, 1, 0, 0, 5)
    assert i2p.inc() is True
    assert i2p.fields == (2, 0, 0, 0, 1, 0, 10)
    assert i2p.inc() is True
    assert i2p.fields == (1, 1, 0, 0, 1, 0, 11)


def test_i2_pow_raise():
    try:
        i2p = I2Pow(-1,1,3)
        assert False
    except ValueError:
        pass

    try:
        i2p = I2Pow(1,-1,3)
        assert False
    except ValueError:
        pass

    try:
        i2p = I2Pow(2,1,3)
        assert False
    except ValueError:
        pass

    try:
        i2p = I2Pow(1,2,3)
        assert False
    except ValueError:
        pass
