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


def test_iter_pow1_inc_l0():
    indexes = np.zeros(3, int)
    result = iter_pow1_inc(indexes)
    assert (indexes == (0,0,0)).all()
    assert result is False


def test_iter_pow1_inc_l1():
    indexes = np.array([1,0,0])
    result = iter_pow1_inc(indexes)
    assert (indexes == (0,1,0)).all()
    assert result is True
    result = iter_pow1_inc(indexes)
    assert (indexes == (0,0,1)).all()
    assert result is True
    result = iter_pow1_inc(indexes)
    assert (indexes == (1,0,0)).all()
    assert result is False


def test_iter_pow1_inc_l2():
    indexes = np.array([2,0,0])
    result = iter_pow1_inc(indexes)
    assert (indexes == (1,1,0)).all()
    assert result is True
    result = iter_pow1_inc(indexes)
    assert (indexes == (1,0,1)).all()
    assert result is True
    result = iter_pow1_inc(indexes)
    assert (indexes == (0,2,0)).all()
    assert result is True
    result = iter_pow1_inc(indexes)
    assert (indexes == (0,1,1)).all()
    assert result is True
    result = iter_pow1_inc(indexes)
    assert (indexes == (0,0,2)).all()
    assert result is True
    result = iter_pow1_inc(indexes)
    assert (indexes == (2,0,0)).all()
    assert result is False


def test_iter_pow2_l0l0():
    ip2 = IterPow2(0, 0, 3)
    assert ip2.fields == (0, 0, 0, 0, 0, 0, 0)
    assert ip2.inc() is False
    assert ip2.fields == (0, 0, 0, 0, 0, 0, 0)


def test_iter_pow2_l1l0():
    ip2 = IterPow2(1, 0, 3)
    assert ip2.fields == (1, 0, 0, 0, 0, 0, 0)
    assert ip2.inc() is True
    assert ip2.fields == (0, 1, 0, 0, 0, 0, 3)
    assert ip2.inc() is True
    assert ip2.fields == (0, 0, 1, 0, 0, 0, 6)
    assert ip2.inc() is False
    assert ip2.fields == (1, 0, 0, 0, 0, 0, 0)


def test_iter_pow2_l2l1():
    ip2 = IterPow2(2, 1, 10)
    assert ip2.fields == (2, 0, 0, 1, 0, 0, 0)
    assert ip2.inc() is True
    assert ip2.fields == (2, 0, 0, 0, 1, 0, 1)
    assert ip2.inc() is True
    assert ip2.fields == (2, 0, 0, 0, 0, 1, 2)
    assert ip2.inc() is True
    assert ip2.fields == (1, 1, 0, 1, 0, 0, 10)
    assert ip2.inc() is True
    assert ip2.fields == (1, 1, 0, 0, 1, 0, 11)
    assert ip2.inc() is True
    assert ip2.fields == (1, 1, 0, 0, 0, 1, 12)
    assert ip2.inc() is True
    assert ip2.fields == (1, 0, 1, 1, 0, 0, 20)
    assert ip2.inc() is True
    assert ip2.fields == (1, 0, 1, 0, 1, 0, 21)


def test_iter_pow2_raise():
    try:
        ip2 = IterPow2(-1,1,3)
        assert False
    except ValueError:
        pass

    try:
        ip2 = IterPow2(1,-1,3)
        assert False
    except ValueError:
        pass

    try:
        ip2 = IterPow2(2,1,3)
        assert False
    except ValueError:
        pass

    try:
        ip2 = IterPow2(1,2,3)
        assert False
    except ValueError:
        pass



def get_test_itergb2():
    # This is a rather random basis specification for two centers.
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 0, 1, 1, 1, 1])
    nprims = np.array([2, 3, 3, 5, 5, 5, 7])
    shell_types = np.array([2, 1, 0, -2, 3, 0, 1])
    alphas = np.random.uniform(-1, 1, nprims.sum())
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())
    gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    return IterGB2(gobasis), gobasis


def test_itergb2_inc_shell():
    i2, gobasis = get_test_itergb2()
    i2.update_shell()
    i2.update_prim()
    assert i2.private_fields == (0, 0,   2, 2, 0, 0, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[0]*gobasis.con_coeffs[0], 2, 2,
        gobasis.alphas[0], gobasis.alphas[0],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        0, 0,
    )
    assert i2.inc_shell() is True
    i2.update_prim()
    assert i2.private_fields == (1, 0,   3, 2, 2, 0, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[0]*gobasis.con_coeffs[2], 1, 2,
        gobasis.alphas[2], gobasis.alphas[0],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        6, 0,
    )
    assert i2.inc_shell() is True
    i2.update_prim()
    assert i2.private_fields == (1, 1,   3, 3, 2, 2, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[2]*gobasis.con_coeffs[2], 1, 1,
        gobasis.alphas[2], gobasis.alphas[2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        6, 6,
    )
    assert i2.inc_shell() is True
    i2.update_prim()
    assert i2.private_fields == (2, 0,   3, 2, 5, 0, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[5]*gobasis.con_coeffs[0], 0, 2,
        gobasis.alphas[5], gobasis.alphas[0],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9, 0,
    )
    assert i2.inc_shell() is True
    i2.update_prim()
    assert i2.private_fields == (2, 1,   3, 3, 5, 2, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[5]*gobasis.con_coeffs[2], 0, 1,
        gobasis.alphas[5], gobasis.alphas[2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9, 6,
    )
    assert i2.inc_shell() is True
    i2.update_prim()
    assert i2.private_fields == (2, 2,   3, 3, 5, 5, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[5]*gobasis.con_coeffs[5], 0, 0,
        gobasis.alphas[5], gobasis.alphas[5],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9, 9,
    )
    assert i2.inc_shell() is True
    i2.update_prim()
    assert i2.private_fields == (3, 0,   5, 2, 8, 0, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[8]*gobasis.con_coeffs[0], -2, 2,
        gobasis.alphas[8], gobasis.alphas[0],
        gobasis.centers[1,0], gobasis.centers[1,1], gobasis.centers[1,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        10, 0,
    )


def test_itergb2_inc_prim():
    i2, gobasis = get_test_itergb2()
    i2.update_shell()
    i2.update_prim()
    assert i2.private_fields == (0, 0,   2, 2, 0, 0, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[0]*gobasis.con_coeffs[0], 2, 2,
        gobasis.alphas[0], gobasis.alphas[0],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        0, 0,
    )
    assert i2.inc_shell() is True
    assert i2.inc_shell() is True
    assert i2.inc_shell() is True
    assert i2.inc_shell() is True
    i2.update_prim()
    assert i2.private_fields == (2, 1,   3, 3, 5, 2, 0, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[5]*gobasis.con_coeffs[2], 0, 1,
        gobasis.alphas[5], gobasis.alphas[2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9, 6,
    )
    assert i2.inc_prim() is True
    assert i2.private_fields == (2, 1,   3, 3, 5, 2, 0, 1)
    assert i2.public_fields == (
        gobasis.con_coeffs[5]*gobasis.con_coeffs[3], 0, 1,
        gobasis.alphas[5], gobasis.alphas[3],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9, 6,
    )
    assert i2.inc_prim() is True
    assert i2.private_fields == (2, 1,   3, 3, 5, 2, 0, 2)
    assert i2.public_fields == (
        gobasis.con_coeffs[5]*gobasis.con_coeffs[4], 0, 1,
        gobasis.alphas[5], gobasis.alphas[4],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9, 6,
    )
    assert i2.inc_prim() is True
    assert i2.private_fields == (2, 1,   3, 3, 5, 2, 1, 0)
    assert i2.public_fields == (
        gobasis.con_coeffs[6]*gobasis.con_coeffs[2], 0, 1,
        gobasis.alphas[6], gobasis.alphas[2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9, 6,
    )


def test_itergb2_store():
    i2, gobasis = get_test_itergb2()
    i2.update_shell()
    i2.update_prim()
    #
    max_shell_nbasis = get_shell_nbasis(gobasis.max_shell_type)
    tmp = np.random.uniform(-1, 1, (max_shell_nbasis, max_shell_nbasis))
    work = tmp + tmp.T
    output = np.zeros((29, 29), float)
    i2.store(work, output)
    assert abs(output[:6,:6] - work[:6,:6]).max() < 1e-10
    assert abs(output[6:,:]).max() == 0
    assert abs(output[:,6:]).max() == 0
    #
    work = np.random.uniform(-1, 1, (max_shell_nbasis, max_shell_nbasis))
    output[:] = 0
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
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
    i2.inc_shell()
    output[:] = 0
    i2.store(work, output)
    assert abs(output[10:15,9:10] - work[:5,:1]).max() < 1e-10
    assert abs(output[9:10,10:15].T - work[:5,:1]).max() < 1e-10
    assert abs(output[:9,:]).max() == 0
    assert abs(output[15:,:]).max() == 0
    assert abs(output[:,:9]).max() == 0
    assert abs(output[:,15:]).max() == 0
    assert abs(output[:10,:10]).max() == 0
    assert abs(output[11:,11:]).max() == 0
