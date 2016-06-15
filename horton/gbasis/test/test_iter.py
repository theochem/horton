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


import numpy as np
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


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


def test_iter_pow1_l0():
    ip1 = IterPow1(0)
    assert ip1.fields == (0, 0, 0, 0)
    assert ip1.inc() is False
    assert ip1.fields == (0, 0, 0, 0)


def test_iter_pow1_l1():
    ip1 = IterPow1(1)
    assert ip1.fields == (1, 0, 0, 0)
    assert ip1.inc() is True
    assert ip1.fields == (0, 1, 0, 1)
    assert ip1.inc() is True
    assert ip1.fields == (0, 0, 1, 2)
    assert ip1.inc() is False
    assert ip1.fields == (1, 0, 0, 0)


def test_iter_pow1_l2():
    ip1 = IterPow1(2)
    assert ip1.fields == (2, 0, 0, 0)
    assert ip1.inc() is True
    assert ip1.fields == (1, 1, 0, 1)
    assert ip1.inc() is True
    assert ip1.fields == (1, 0, 1, 2)
    assert ip1.inc() is True
    assert ip1.fields == (0, 2, 0, 3)
    assert ip1.inc() is True
    assert ip1.fields == (0, 1, 1, 4)
    assert ip1.inc() is True
    assert ip1.fields == (0, 0, 2, 5)
    assert ip1.inc() is False
    assert ip1.fields == (2, 0, 0, 0)


def test_iter_pow2_l0l0():
    ip2 = IterPow2(0, 0)
    assert ip2.fields == (0, 0, 0, 0, 0, 0, 0, 0, 0)
    assert ip2.inc() is False
    assert ip2.fields == (0, 0, 0, 0, 0, 0, 0, 0, 0)


def test_iter_pow2_l1l0():
    ip2 = IterPow2(1, 0)
    assert ip2.fields == (1, 0, 0, 0, 0, 0, 0, 0, 0)
    assert ip2.inc() is True
    assert ip2.fields == (0, 1, 0, 0, 0, 0, 1, 1, 0)
    assert ip2.inc() is True
    assert ip2.fields == (0, 0, 1, 0, 0, 0, 2, 2, 0)
    assert ip2.inc() is False
    assert ip2.fields == (1, 0, 0, 0, 0, 0, 0, 0, 0)


def test_iter_pow2_l2l1():
    ip2 = IterPow2(2, 1)
    assert ip2.fields == (2, 0, 0, 1, 0, 0, 0, 0, 0)
    assert ip2.inc() is True
    assert ip2.fields == (2, 0, 0, 0, 1, 0, 1, 0, 1)
    assert ip2.inc() is True
    assert ip2.fields == (2, 0, 0, 0, 0, 1, 2, 0, 2)
    assert ip2.inc() is True
    assert ip2.fields == (1, 1, 0, 1, 0, 0, 3, 1, 0)
    assert ip2.inc() is True
    assert ip2.fields == (1, 1, 0, 0, 1, 0, 4, 1, 1)
    assert ip2.inc() is True
    assert ip2.fields == (1, 1, 0, 0, 0, 1, 5, 1, 2)
    assert ip2.inc() is True
    assert ip2.fields == (1, 0, 1, 1, 0, 0, 6, 2, 0)
    assert ip2.inc() is True
    assert ip2.fields == (1, 0, 1, 0, 1, 0, 7, 2, 1)


def test_iter_pow2_raise():
    with assert_raises(ValueError):
        ip2 = IterPow2(-1,1)
    with assert_raises(ValueError):
        ip2 = IterPow2(1,-1)


def get_itergb1():
    # This is a rather random basis specification for two centers.
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 0, 1, 1, 1, 1])
    nprims = np.array([2, 3, 3, 5, 5, 5, 7])
    shell_types = np.array([2, 1, 0, -2, 3, 0, 1])
    alphas = np.random.uniform(0.5, 2, nprims.sum())
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())
    gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    return IterGB1(gobasis), gobasis


def test_itergb1_inc_shell():
    i1, gobasis = get_itergb1()
    i1.update_shell()
    i1.update_prim()
    assert i1.private_fields == (0,   2, 0, 0)
    assert i1.public_fields == (
        gobasis.con_coeffs[0], 2,
        gobasis.alphas[0],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        0,
    )
    assert i1.inc_shell() is True
    i1.update_prim()
    assert i1.private_fields == (1,   3, 2, 0)
    assert i1.public_fields == (
        gobasis.con_coeffs[2], 1,
        gobasis.alphas[2],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        6,
    )
    assert i1.inc_shell() is True
    i1.update_prim()
    assert i1.private_fields == (2,   3, 5, 0)
    assert i1.public_fields == (
        gobasis.con_coeffs[5], 0,
        gobasis.alphas[5],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9,
    )
    assert i1.inc_shell() is True
    i1.update_prim()
    assert i1.private_fields == (3,   5, 8, 0)
    assert i1.public_fields == (
        gobasis.con_coeffs[8], -2,
        gobasis.alphas[8],
        gobasis.centers[1,0], gobasis.centers[1,1], gobasis.centers[1,2],
        10,
    )
    assert i1.inc_shell() is True
    assert i1.inc_shell() is True
    assert i1.inc_shell() is True
    assert i1.inc_shell() is False



def test_itergb1_inc_prim():
    i1, gobasis = get_itergb1()
    i1.update_shell()
    i1.update_prim()
    assert i1.private_fields == (0,   2, 0, 0)
    assert i1.public_fields == (
        gobasis.con_coeffs[0], 2,
        gobasis.alphas[0],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        0,
    )
    assert i1.inc_shell() is True
    assert i1.inc_shell() is True
    i1.update_prim()
    assert i1.private_fields == (2,   3, 5, 0)
    assert i1.public_fields == (
        gobasis.con_coeffs[5], 0,
        gobasis.alphas[5],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9
    )
    assert i1.inc_prim() is True
    assert i1.private_fields == (2,  3, 5, 1)
    assert i1.public_fields == (
        gobasis.con_coeffs[6], 0,
        gobasis.alphas[6],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9,
    )
    assert i1.inc_prim() is True
    assert i1.private_fields == (2,  3, 5, 2)
    assert i1.public_fields == (
        gobasis.con_coeffs[7], 0,
        gobasis.alphas[7],
        gobasis.centers[0,0], gobasis.centers[0,1], gobasis.centers[0,2],
        9,
    )
    assert i1.inc_prim() is False


def test_itergb1_store():
    i1, gobasis = get_itergb1()
    i1.update_shell()
    i1.update_prim()
    #
    max_shell_nbasis = get_shell_nbasis(gobasis.max_shell_type)
    work = np.random.uniform(-1, 1, 6)
    output = np.zeros(29, float)
    i1.store(work, output)
    assert abs(output[:6] - work).max() < 1e-10
    #
    work = np.random.uniform(-1, 1, 5)
    output[:] = 0
    i1.inc_shell()
    i1.inc_shell()
    i1.inc_shell()
    i1.store(work, output)
    assert abs(output[10:15] - work).max() < 1e-10
    assert abs(output[:10]).max() == 0
    assert abs(output[15:]).max() == 0



def get_itergb2():
    # This is a rather random basis specification for two centers.
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 0, 1, 1, 1, 1])
    nprims = np.array([2, 3, 3, 5, 5, 5, 7])
    shell_types = np.array([2, 1, 0, -2, 3, 0, 1])
    alphas = np.random.uniform(0.5, 2, nprims.sum())
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())
    gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    return IterGB2(gobasis), gobasis


def test_itergb2_inc_shell():
    i2, gobasis = get_itergb2()
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
    i2, gobasis = get_itergb2()
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
    i2, gobasis = get_itergb2()
    i2.update_shell()
    i2.update_prim()
    #
    max_shell_nbasis = get_shell_nbasis(gobasis.max_shell_type)
    tmp = np.random.uniform(-1, 1, (6, 6))
    work = tmp + tmp.T
    output = np.zeros((29, 29), float)
    i2.store(work, output)
    assert abs(output[:6,:6] - work).max() < 1e-10
    assert abs(output[6:,:]).max() == 0
    assert abs(output[:,6:]).max() == 0
    #
    work = np.random.uniform(-1, 1, (5, 3))
    output[:] = 0
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.inc_shell()
    i2.store(work, output)
    assert abs(output[10:15,6:9] - work).max() < 1e-10
    assert abs(output[6:9,10:15].T - work).max() < 1e-10
    assert abs(output[:6,:]).max() == 0
    assert abs(output[15:,:]).max() == 0
    assert abs(output[:,:6]).max() == 0
    assert abs(output[:,15:]).max() == 0
    assert abs(output[:10,:10]).max() == 0
    assert abs(output[9:,9:]).max() == 0
    #
    i2.inc_shell()
    output[:] = 0
    work = np.random.uniform(-1, 1, (5, 1))
    i2.store(work, output)
    assert abs(output[10:15,9:10] - work).max() < 1e-10
    assert abs(output[9:10,10:15].T - work).max() < 1e-10
    assert abs(output[:9,:]).max() == 0
    assert abs(output[15:,:]).max() == 0
    assert abs(output[:,:9]).max() == 0
    assert abs(output[:,15:]).max() == 0
    assert abs(output[:10,:10]).max() == 0
    assert abs(output[11:,11:]).max() == 0


def test_itergb4_idea():
    # # This is just a test for a loop structure idea.
    # # Note: physicists notation is used.
    # a counter for each element of a four-index operator with 5 DOF:
    tboc = np.zeros((5, 5, 5, 5), int)
    # quadruple loop over all unique elements:
    for i in xrange(5):
        for j in xrange(i+1): # (i>=j)
            for k in xrange(i+1): # (i>=k)
                if i==j:
                    # this is a special case: when i==j, the permutation of the
                    # last two indexes is also a symmetry. Hence we have (k>=l)
                    lm = k
                else:
                    # this is the regular case: (j>=l)
                    lm = j
                for l in xrange(lm+1):
                    # add 1 to tboc element and all relevant symmetries
                    tboc[i,j,k,l] += 1
                    tboc[j,i,l,k] += 1
                    tboc[k,l,i,j] += 1
                    tboc[l,k,j,i] += 1
                    tboc[k,j,i,l] += 1
                    tboc[l,i,j,k] += 1
                    tboc[i,l,k,j] += 1
                    tboc[j,k,l,i] += 1
    # The counter should be equal to the number of relevant permutations that
    # leave the indices unchanged.
    for i in xrange(5):
        for j in xrange(5):
            for k in xrange(5):
                for l in xrange(5):
                    expected = 1+sum([
                        (i,j,k,l) == (j,i,l,k),
                        (i,j,k,l) == (k,l,i,j),
                        (i,j,k,l) == (l,k,j,i),
                        (i,j,k,l) == (i,l,k,j),
                        (i,j,k,l) == (j,k,l,i),
                        (i,j,k,l) == (k,j,i,l),
                        (i,j,k,l) == (l,i,j,k),
                    ])
                    good = tboc[i,j,k,l] == expected
                    if not good:
                        raise AssertionError('(%i,%i,%i,%i): %i != %i' % (i,j,k,l,tboc[i,j,k,l],expected))


def get_itergb4():
    # This is a rather random basis specification for two centers.
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 1, 1, 1])
    nprims = np.array([2, 3, 5, 2, 7])
    shell_types = np.array([2, 1, -2, 3, 0])
    alphas = np.random.uniform(0.5, 2, nprims.sum())
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())
    gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    return IterGB4(gobasis), gobasis


def test_itergb4_idea_inc_shell():
    # # This is a test for the shell loop structure of IterGB4.
    # # Note: physicists notation is used.
    i4, gobasis = get_itergb4()
    N = gobasis.nshell
    # a counter for each element of a four-index operator with 5 DOF:
    tboc = np.zeros((N, N, N, N), int)
    # quadruple loop over all unique elements:
    i4.update_shell()
    while True:
        i, j, k, l = i4.private_fields[:4]
        tboc[i,j,k,l] += 1
        tboc[j,i,l,k] += 1
        tboc[k,l,i,j] += 1
        tboc[l,k,j,i] += 1
        tboc[k,j,i,l] += 1
        tboc[l,i,j,k] += 1
        tboc[i,l,k,j] += 1
        tboc[j,k,l,i] += 1
        if not i4.inc_shell():
            break
    # The counter should be equal to the number of relevant permutations that
    # leave the indices unchanged.
    for i in xrange(N):
        for j in xrange(N):
            for k in xrange(N):
                for l in xrange(N):
                    expected = 1+sum([
                        (i,j,k,l) == (j,i,l,k),
                        (i,j,k,l) == (k,l,i,j),
                        (i,j,k,l) == (l,k,j,i),
                        (i,j,k,l) == (i,l,k,j),
                        (i,j,k,l) == (j,k,l,i),
                        (i,j,k,l) == (k,j,i,l),
                        (i,j,k,l) == (l,i,j,k),
                    ])
                    good = tboc[i,j,k,l] == expected
                    if not good:
                        raise AssertionError('(%i,%i,%i,%i): %i != %i' % (i,j,k,l,tboc[i,j,k,l],expected))


def test_itergb4_inc_shell():
    i4, gobasis = get_itergb4()
    oprims = np.zeros(gobasis.nshell, int)
    basis_offsets = np.zeros(gobasis.nshell, int)
    for i in xrange(1, gobasis.nshell):
        oprims[i] = oprims[i-1] + gobasis.nprims[i-1]
        basis_offsets[i] = basis_offsets[i-1] + get_shell_nbasis(gobasis.shell_types[i-1])

    def check_fields(is0, is1, is2, is3):
        assert i4.private_fields == (
            is0, is1, is2, is3,
            gobasis.nprims[is0], gobasis.nprims[is1], gobasis.nprims[is2], gobasis.nprims[is3],
            oprims[is0], oprims[is1], oprims[is2], oprims[is3],
            0, 0, 0, 0
        )
        f1 = np.array(list(i4.public_fields))
        f2 = np.array([
            gobasis.con_coeffs[oprims[is0]]*gobasis.con_coeffs[oprims[is1]]*gobasis.con_coeffs[oprims[is2]]*gobasis.con_coeffs[oprims[is3]],
            gobasis.shell_types[is0], gobasis.shell_types[is1], gobasis.shell_types[is2], gobasis.shell_types[is3],
            gobasis.alphas[oprims[is0]], gobasis.alphas[oprims[is1]], gobasis.alphas[oprims[is2]], gobasis.alphas[oprims[is3]],
            gobasis.centers[gobasis.shell_map[is0],0], gobasis.centers[gobasis.shell_map[is0],1], gobasis.centers[gobasis.shell_map[is0],2],
            gobasis.centers[gobasis.shell_map[is1],0], gobasis.centers[gobasis.shell_map[is1],1], gobasis.centers[gobasis.shell_map[is1],2],
            gobasis.centers[gobasis.shell_map[is2],0], gobasis.centers[gobasis.shell_map[is2],1], gobasis.centers[gobasis.shell_map[is2],2],
            gobasis.centers[gobasis.shell_map[is3],0], gobasis.centers[gobasis.shell_map[is3],1], gobasis.centers[gobasis.shell_map[is3],2],
            basis_offsets[is0], basis_offsets[is1], basis_offsets[is2], basis_offsets[is3]
        ])
        assert abs(f1 - f2).max() < 1e-10

    i4.update_shell()
    i4.update_prim()
    check_fields(0, 0, 0, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(1, 0, 0, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(1, 0, 1, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(1, 1, 0, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(1, 1, 1, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(1, 1, 1, 1)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 0, 0, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 0, 1, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 0, 2, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 1, 0, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 1, 0, 1)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 1, 1, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 1, 1, 1)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 1, 2, 0)

    assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 1, 2, 1)


def test_itergb4_inc_prim():
    i4, gobasis = get_itergb4()
    oprims = np.zeros(gobasis.nshell, int)
    basis_offsets = np.zeros(gobasis.nshell, int)
    for i in xrange(1, gobasis.nshell):
        oprims[i] = oprims[i-1] + gobasis.nprims[i-1]
        basis_offsets[i] = basis_offsets[i-1] + get_shell_nbasis(gobasis.shell_types[i-1])

    def check_fields(is0, is1, is2, is3, ip0, ip1, ip2, ip3):
        assert i4.private_fields == (
            is0, is1, is2, is3,
            gobasis.nprims[is0], gobasis.nprims[is1], gobasis.nprims[is2], gobasis.nprims[is3],
            oprims[is0], oprims[is1], oprims[is2], oprims[is3],
            ip0, ip1, ip2, ip3
        )
        f1 = np.array(list(i4.public_fields))
        f2 = np.array([
            gobasis.con_coeffs[oprims[is0]+ip0]*gobasis.con_coeffs[oprims[is1]+ip1]*gobasis.con_coeffs[oprims[is2]+ip2]*gobasis.con_coeffs[oprims[is3]+ip3],
            gobasis.shell_types[is0], gobasis.shell_types[is1], gobasis.shell_types[is2], gobasis.shell_types[is3],
            gobasis.alphas[oprims[is0]+ip0], gobasis.alphas[oprims[is1]+ip1], gobasis.alphas[oprims[is2]+ip2], gobasis.alphas[oprims[is3]+ip3],
            gobasis.centers[gobasis.shell_map[is0],0], gobasis.centers[gobasis.shell_map[is0],1], gobasis.centers[gobasis.shell_map[is0],2],
            gobasis.centers[gobasis.shell_map[is1],0], gobasis.centers[gobasis.shell_map[is1],1], gobasis.centers[gobasis.shell_map[is1],2],
            gobasis.centers[gobasis.shell_map[is2],0], gobasis.centers[gobasis.shell_map[is2],1], gobasis.centers[gobasis.shell_map[is2],2],
            gobasis.centers[gobasis.shell_map[is3],0], gobasis.centers[gobasis.shell_map[is3],1], gobasis.centers[gobasis.shell_map[is3],2],
            basis_offsets[is0], basis_offsets[is1], basis_offsets[is2], basis_offsets[is3]
        ])
        assert abs(f1 - f2).max() < 1e-10

    i4.update_shell()
    i4.update_prim()
    check_fields(0, 0, 0, 0,  0, 0, 0, 0)

    for i in xrange(10):
        assert i4.inc_shell() is True
    i4.update_prim()
    check_fields(2, 1, 0, 1,  0, 0, 0, 0)

    assert i4.inc_prim() is True
    check_fields(2, 1, 0, 1,  0, 0, 0, 1)

    assert i4.inc_prim() is True
    check_fields(2, 1, 0, 1,  0, 0, 0, 2)

    assert i4.inc_prim() is True
    check_fields(2, 1, 0, 1,  0, 0, 1, 0)

    for i in xrange(5*3*2*3-4):
        assert i4.inc_prim() is True

    assert i4.inc_prim() is False
    check_fields(2, 1, 0, 1,  0, 0, 0, 0)


def test_itergb4_store():
    i4, gobasis = get_itergb4()
    i4.update_shell()
    i4.update_prim()
    #
    max_shell_nbasis = get_shell_nbasis(gobasis.max_shell_type)
    work = np.random.uniform(-1, 1, (6,)*4)
    output = np.zeros((25, 25, 25, 25), float)
    #
    i4.store(work, output)
    # it is hard to predict how an asymmetric work array will be stored when
    # the indexes are (partly) on the diagonal
    assert (output[:6,:6,:6,:6] != 0).any()
    assert abs(output[6:,:,:,:]).max() == 0
    assert abs(output[:,6:,:,:]).max() == 0
    assert abs(output[:,:,6:,:]).max() == 0
    assert abs(output[:,:,:,6:]).max() == 0
    #
    while i4.private_fields[:4] != (3, 2, 1, 0):
        assert i4.inc_shell() is True
    work = np.random.uniform(-1, 1, (10, 5, 3, 6))
    output[:] = 0
    mask = np.ones(output.shape, bool)
    i4.store(work, output)
    slices = [slice(14,24), slice(9,14), slice(6,9), slice(0,6)]
    permutations = [
        (0, 1, 2, 3), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0),
        (0, 3, 2, 1), (1, 2, 3, 0), (2, 1, 0, 3), (3, 0, 1, 2),
    ]
    for i0, i1, i2, i3 in permutations:
        assert abs(output[slices[i0],slices[i1],slices[i2],slices[i3]] - work.transpose(i0,i1,i2,i3)).max() < 1e-10
        mask[slices[i0],slices[i1],slices[i2],slices[i3]] = False
    assert abs(output[mask]).max() == 0.0
