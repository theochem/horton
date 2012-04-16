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
    # This is a rather random basis specification for two centers.
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 1, 1])
    ncons = np.array([1, 2, 3, 1])
    nexps = np.array([2, 3, 5, 7])
    con_types = np.array([2, 1, 0, -2, 3, 0, 1])
    exponents = np.random.uniform(-1, 1, nexps.sum())
    con_coeffs = np.random.uniform(-1, 1, np.dot(nexps, ncons))
    nbasis = sum(get_con_nbasis(con_type) for con_type in con_types)
    assert nbasis == 29
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
    i2.update_shell()
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (0, 0,   1, 1, 0, 0, 0, 0,   2, 2, 0, 0, 0, 0,   0, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[0], 2, 2,
        i2.exponents[0], i2.exponents[0],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        0, 0,
    )
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (1, 0,   2, 1, 1, 0, 0, 0,   3, 2, 2, 0, 0, 0,   2, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[2], 1, 2,
        i2.exponents[2], i2.exponents[0],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        6, 0,
    )
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (1, 1,   2, 2, 1, 1, 0, 0,   3, 3, 2, 2, 0, 0,   2, 2)
    assert i2.public_fields == (
        i2.con_coeffs[2]*i2.con_coeffs[2], 1, 1,
        i2.exponents[2], i2.exponents[2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        6, 6,
    )
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (2, 0,   3, 1, 3, 0, 0, 0,   5, 2, 5, 0, 0, 0,   8, 0)
    assert i2.public_fields == (
        i2.con_coeffs[8]*i2.con_coeffs[0], -2, 2,
        i2.exponents[5], i2.exponents[0],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        10, 0,
    )
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (2, 1,   3, 2, 3, 1, 0, 0,   5, 3, 5, 2, 0, 0,   8, 2)
    assert i2.public_fields == (
        i2.con_coeffs[2]*i2.con_coeffs[8], -2, 1,
        i2.exponents[5], i2.exponents[2],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        10, 6,
    )
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (2, 2,   3, 3, 3, 3, 0, 0,   5, 5, 5, 5, 0, 0,   8, 8)
    assert i2.public_fields == (
        i2.con_coeffs[8]*i2.con_coeffs[8], -2, -2,
        i2.exponents[5], i2.exponents[5],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        10, 10,
    )
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (3, 0,   1, 1, 6, 0, 0, 0,   7, 2, 10, 0, 0, 0,   23, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[23], 1, 2,
        i2.exponents[10], i2.exponents[0],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        26, 0,
    )


def test_i2gob_inc_con():
    i2 = get_test_i2gob()[0]
    i2.update_shell()
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (0, 0,   1, 1, 0, 0, 0, 0,   2, 2, 0, 0, 0, 0,   0, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[0], 2, 2,
        i2.exponents[0], i2.exponents[0],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        0, 0,
    )
    assert i2.inc_con() is False
    i2.update_exp()
    assert i2.private_fields == (0, 0,   1, 1, 0, 0, 0, 0,   2, 2, 0, 0, 0, 0,   0, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[0], 2, 2,
        i2.exponents[0], i2.exponents[0],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        0, 0,
    )
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (1, 0,   2, 1, 1, 0, 0, 0,   3, 2, 2, 0, 0, 0,   2, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[2], 1, 2,
        i2.exponents[2], i2.exponents[0],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        6, 0,
    )
    assert i2.inc_con() is True
    i2.update_exp()
    assert i2.private_fields == (1, 0,   2, 1, 1, 0, 1, 0,   3, 2, 2, 0, 0, 0,   2, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[3], 0, 2,
        i2.exponents[2], i2.exponents[0],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        9, 0,
    )
    assert i2.inc_con() is False
    i2.update_exp()
    assert i2.private_fields == (1, 0,   2, 1, 1, 0, 0, 0,   3, 2, 2, 0, 0, 0,   2, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[2], 1, 2,
        i2.exponents[2], i2.exponents[0],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        6, 0,
    )
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (1, 1,   2, 2, 1, 1, 0, 0,   3, 3, 2, 2, 0, 0,   2, 2)
    assert i2.public_fields == (
        i2.con_coeffs[2]*i2.con_coeffs[2], 1, 1,
        i2.exponents[2], i2.exponents[2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        6, 6,
    )
    assert i2.inc_con() is True
    i2.update_exp()
    assert i2.private_fields == (1, 1,   2, 2, 1, 1, 0, 1,   3, 3, 2, 2, 0, 0,   2, 2)
    assert i2.public_fields == (
        i2.con_coeffs[3]*i2.con_coeffs[2], 1, 0,
        i2.exponents[2], i2.exponents[2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        6, 9,
    )
    assert i2.inc_con() is True
    i2.update_exp()
    assert i2.private_fields == (1, 1,   2, 2, 1, 1, 1, 0,   3, 3, 2, 2, 0, 0,   2, 2)
    assert i2.public_fields == (
        i2.con_coeffs[2]*i2.con_coeffs[3], 0, 1,
        i2.exponents[2], i2.exponents[2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        9, 6,
    )
    assert i2.inc_con() is True
    i2.update_exp()
    assert i2.private_fields == (1, 1,   2, 2, 1, 1, 1, 1,   3, 3, 2, 2, 0, 0,   2, 2)
    assert i2.public_fields == (
        i2.con_coeffs[3]*i2.con_coeffs[3], 0, 0,
        i2.exponents[2], i2.exponents[2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        9, 9,
    )


def test_i2gob_inc_exp():
    i2 = get_test_i2gob()[0]
    i2.update_shell()
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (0, 0,   1, 1, 0, 0, 0, 0,   2, 2, 0, 0, 0, 0,   0, 0)
    assert i2.public_fields == (
        i2.con_coeffs[0]*i2.con_coeffs[0], 2, 2,
        i2.exponents[0], i2.exponents[0],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        0, 0,
    )
    assert i2.inc_shell() is True
    assert i2.inc_shell() is True
    assert i2.inc_shell() is True
    assert i2.inc_shell() is True
    i2.update_con()
    i2.update_exp()
    assert i2.private_fields == (2, 1,   3, 2, 3, 1, 0, 0,   5, 3, 5, 2, 0, 0,   8, 2)
    assert i2.public_fields == (
        i2.con_coeffs[2]*i2.con_coeffs[8], -2, 1,
        i2.exponents[5], i2.exponents[2],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        10, 6,
    )
    assert i2.inc_exp() is True
    assert i2.private_fields == (2, 1,   3, 2, 3, 1, 0, 0,   5, 3, 5, 2, 0, 1,   8, 2)
    assert i2.public_fields == (
        i2.con_coeffs[4]*i2.con_coeffs[8], -2, 1,
        i2.exponents[5], i2.exponents[3],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        10, 6,
    )
    assert i2.inc_exp() is True
    assert i2.private_fields == (2, 1,   3, 2, 3, 1, 0, 0,   5, 3, 5, 2, 0, 2,   8, 2)
    assert i2.public_fields == (
        i2.con_coeffs[6]*i2.con_coeffs[8], -2, 1,
        i2.exponents[5], i2.exponents[4],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        10, 6,
    )
    assert i2.inc_exp() is True
    assert i2.private_fields == (2, 1,   3, 2, 3, 1, 0, 0,   5, 3, 5, 2, 1, 0,   8, 2)
    assert i2.public_fields == (
        i2.con_coeffs[2]*i2.con_coeffs[11], -2, 1,
        i2.exponents[6], i2.exponents[2],
        i2.centers[1,0], i2.centers[1,1], i2.centers[1,2],
        i2.centers[0,0], i2.centers[0,1], i2.centers[0,2],
        10, 6,
    )


def test_i2gob_store():
    i2 = get_test_i2gob()[0]
    i2.update_shell()
    i2.update_con()
    i2.update_exp()
    #
    tmp = np.random.uniform(-1, 1, (i2.max_nbasis, i2.max_nbasis))
    work = tmp + tmp.T
    output = np.zeros((29, 29), float)
    i2.store(work, output)
    assert abs(output[:6,:6] - work[:6,:6]).max() < 1e-10
    assert abs(output[6:,:]).max() == 0
    assert abs(output[:,6:]).max() == 0
    #
    work = np.random.uniform(-1, 1, (i2.max_nbasis, i2.max_nbasis))
    output[:] = 0
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.update_con()
    i2.store(work, output)
    assert abs(output[10:15,6:9] - work[:5,:3]).max() < 1e-10
    assert abs(output[6:9,10:15].T - work[:5,:3]).max() < 1e-10
    assert abs(output[:6,:]).max() == 0
    assert abs(output[15:,:]).max() == 0
    assert abs(output[:,:6]).max() == 0
    assert abs(output[:,15:]).max() == 0
    assert abs(output[:10,:10]).max() == 0
    assert abs(output[9:,9:]).max() == 0
    #
    i2.inc_con()
    output[:] = 0
    i2.store(work, output)
    print (output != 0).astype(int)
    assert abs(output[10:15,9:10] - work[:5,:1]).max() < 1e-10
    assert abs(output[9:10,10:15].T - work[:5,:1]).max() < 1e-10
    assert abs(output[:9,:]).max() == 0
    assert abs(output[15:,:]).max() == 0
    assert abs(output[:,:9]).max() == 0
    assert abs(output[:,15:]).max() == 0
    assert abs(output[:10,:10]).max() == 0
    assert abs(output[11:,11:]).max() == 0


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
    assert i2p.fields == (0, 1, 0, 0, 0, 0, 3)
    assert i2p.inc() is True
    assert i2p.fields == (0, 0, 1, 0, 0, 0, 6)
    assert i2p.inc() is False
    assert i2p.fields == (1, 0, 0, 0, 0, 0, 0)


def test_i2_pow_l2l1():
    i2p = I2Pow(2, 1, 10)
    assert i2p.fields == (2, 0, 0, 1, 0, 0, 0)
    assert i2p.inc() is True
    assert i2p.fields == (2, 0, 0, 0, 1, 0, 1)
    assert i2p.inc() is True
    assert i2p.fields == (2, 0, 0, 0, 0, 1, 2)
    assert i2p.inc() is True
    assert i2p.fields == (1, 1, 0, 1, 0, 0, 10)
    assert i2p.inc() is True
    assert i2p.fields == (1, 1, 0, 0, 1, 0, 11)
    assert i2p.inc() is True
    assert i2p.fields == (1, 1, 0, 0, 0, 1, 12)
    assert i2p.inc() is True
    assert i2p.fields == (1, 0, 1, 1, 0, 0, 20)
    assert i2p.inc() is True
    assert i2p.fields == (1, 0, 1, 0, 1, 0, 21)


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
