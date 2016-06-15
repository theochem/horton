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
from nose.plugins.attrib import attr

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


tfs = {
    2: np.array([
        [-0.5, 0, 0, -0.5, 0, 1.0],
        [0, 0, 1.0, 0, 0, 0],
        [0, 0, 0, 0, 1.0, 0],
        [0.86602540378443864676, 0, 0, -0.86602540378443864676, 0, 0],
        [0, 1.0, 0, 0, 0, 0],
    ]),
    3: np.array([
        [0, 0, -0.67082039324993690892, 0, 0, 0, 0, -0.67082039324993690892, 0, 1.0],
        [-0.61237243569579452455, 0, 0, -0.27386127875258305673, 0, 1.0954451150103322269, 0, 0, 0, 0],
        [0, -0.27386127875258305673, 0, 0, 0, 0, -0.61237243569579452455, 0, 1.0954451150103322269, 0],
        [0, 0, 0.86602540378443864676, 0, 0, 0, 0, -0.86602540378443864676, 0, 0],
        [0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0],
        [0.790569415042094833, 0, 0, -1.0606601717798212866, 0, 0, 0, 0, 0, 0],
        [0, 1.0606601717798212866, 0, 0, 0, 0, -0.790569415042094833, 0, 0, 0],
    ]),
    4: np.array([
        [0.375, 0, 0, 0.21957751641341996535, 0, -0.87831006565367986142, 0, 0, 0, 0, 0.375, 0, -0.87831006565367986142, 0, 1.0],
        [0, 0, -0.89642145700079522998, 0, 0, 0, 0, -0.40089186286863657703, 0, 1.19522860933439364, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -0.40089186286863657703, 0, 0, 0, 0, 0, 0, -0.89642145700079522998, 0, 1.19522860933439364, 0],
        [-0.5590169943749474241, 0, 0, 0, 0, 0.9819805060619657157, 0, 0, 0, 0, 0.5590169943749474241, 0, -0.9819805060619657157, 0, 0],
        [0, -0.42257712736425828875, 0, 0, 0, 0, -0.42257712736425828875, 0, 1.1338934190276816816, 0, 0, 0, 0, 0, 0],
        [0, 0, 0.790569415042094833, 0, 0, 0, 0, -1.0606601717798212866, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1.0606601717798212866, 0, 0, 0, 0, 0, 0, -0.790569415042094833, 0, 0, 0],
        [0.73950997288745200532, 0, 0, -1.2990381056766579701, 0, 0, 0, 0, 0, 0, 0.73950997288745200532, 0, 0, 0, 0],
        [0, 1.1180339887498948482, 0, 0, 0, 0, -1.1180339887498948482, 0, 0, 0, 0, 0, 0, 0, 0],
    ]),
}

def test_cart_pure_s():
    work_cart = np.random.normal(0, 1, 1)
    work_pure = np.random.normal(0, 1, 1)
    cart_to_pure_low(work_cart, work_pure, shell_type=0, nant=1, npost=1)
    assert abs(work_cart - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, 10)
    work_pure = np.random.normal(0, 1, 10)
    cart_to_pure_low(work_cart, work_pure, shell_type=0, nant=2, npost=5)
    assert abs(work_cart - work_pure).max() < 1e-10


def test_cart_pure_p():
    work_cart = np.random.normal(0, 1, 3)
    work_pure = np.random.normal(0, 1, 3)
    cart_to_pure_low(work_cart, work_pure, shell_type=1, nant=1, npost=1)
    assert abs(work_cart[[2,0,1]] - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, (10, 3))
    work_pure = np.random.normal(0, 1, (10, 3))
    cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=1, nant=10, npost=1)
    assert abs(work_cart[:,[2,0,1]] - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, (10, 3, 2))
    work_pure = np.random.normal(0, 1, (10, 3, 2))
    cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=1, nant=10, npost=2)
    assert abs(work_cart[:,[2,0,1],:] - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, (3, 6))
    work_pure = np.random.normal(0, 1, (3, 6))
    cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=1, nant=1, npost=6)
    assert abs(work_cart[[2,0,1],:] - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, (5, 2, 3, 2, 3))
    work_pure = np.random.normal(0, 1, (5, 2, 3, 2, 3))
    cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=1, nant=10, npost=6)
    assert abs(work_cart[:,:,[2,0,1],:,:] - work_pure).max() < 1e-10



def test_cart_pure_d():
    tf = tfs[2]

    work_cart = np.random.normal(0, 1, 6)
    work_pure = np.random.normal(0, 1, 5)
    cart_to_pure_low(work_cart, work_pure, shell_type=2, nant=1, npost=1)
    assert abs(np.dot(tf, work_cart) - work_pure[:5]).max() < 1e-10

    work_cart = np.random.normal(0, 1, (10, 6))
    work_pure = np.random.normal(0, 1, (10, 5))
    cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=2, nant=10, npost=1)
    assert abs(np.dot(work_cart, tf.T) - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, (6, 10))
    work_pure = np.random.normal(0, 1, (5, 10))
    cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=2, nant=1, npost=10)
    assert abs(np.dot(tf, work_cart) - work_pure).max() < 1e-10


def test_cart_pure_g():
    tf = tfs[4]

    work_cart = np.random.normal(0, 1, 15)
    work_pure = np.random.normal(0, 1, 9)
    cart_to_pure_low(work_cart, work_pure, shell_type=4, nant=1, npost=1)
    assert abs(np.dot(tf, work_cart) - work_pure).max() < 1e-10

    work_cart = np.random.normal(0, 1, (3,15))
    work_pure = np.random.normal(0, 1, (3,9))
    cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=4, nant=3, npost=1)
    assert abs(np.dot(work_cart[:,:15], tf.T) - work_pure).max() < 1e-10


def test_gb2_overlap_integral_class():
    max_shell_type = 4
    max_nbasis = get_shell_nbasis(max_shell_type)
    r0 = np.array([2.645617, 0.377945, -0.188973])
    r1 = np.array([1.254878, 0.123456,  0.188973])
    scales0 = np.ones(15, float)
    scales1 = np.ones(10, float)

    gb2oi = GB2OverlapIntegral(max_shell_type)
    assert gb2oi.max_nbasis == max_nbasis

    gb2oi.reset(-4, -3, r0, r1)
    gb2oi.add(1.0, 5.398, 0.320, scales0, scales1)
    gb2oi.add(0.5, 0.123, 0.210, scales0, scales1)
    gb2oi.add(0.7, 1.234, 2.333, scales0, scales1)
    gb2oi.add(0.3, 0.500, 0.500, scales0, scales1)
    work0 = gb2oi.get_work(15,10)
    gb2oi.cart_to_pure()
    work1 = gb2oi.get_work(9,7)
    step0 = np.dot(work0, tfs[3].T)
    step1 = np.dot(tfs[4], step0)
    assert abs(work1 - step1).max() < 1e-10


def test_cart_pure_domain():
    work_cart = np.random.normal(0, 1, (3,70))
    work_pure = np.random.normal(0, 1, (3,70))
    with assert_raises(ValueError):
        cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=get_max_shell_type()+1, nant=1, npost=1)
    with assert_raises(ValueError):
        cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=-1, nant=1, npost=1)
    with assert_raises(ValueError):
        cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=3, nant=0, npost=1)
    with assert_raises(ValueError):
        cart_to_pure_low(work_cart.reshape(-1), work_pure.reshape(-1), shell_type=3, nant=1, npost=0)


@attr('slow')
def test_cart_pure_water_ccpvdz_hf():
    fn_fchk_pure = context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk')
    fn_log_pure = fn_fchk_pure[:-5] + '.log'
    fn_fchk_cart = context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk')
    fn_log_cart = fn_fchk_cart[:-5] + '.log'
    # Also load fchk file to get reordering of matrix elements.
    mol_pure = IOData.from_file(fn_fchk_pure, fn_log_pure)
    mol_cart = IOData.from_file(fn_fchk_cart, fn_log_cart)
    for key in 'olp', 'kin', 'na':
        block_pure = getattr(mol_pure, key)._array[9:14,9:14]
        block_cart = getattr(mol_cart, key)._array[9:15,9:15]
        check_pure = np.dot(np.dot(tfs[2], block_cart), tfs[2].T)
        error = abs(block_pure - check_pure).max()
        assert error < 2e-5

        block_pure = getattr(mol_pure, key)._array[0,9:14]
        block_cart = getattr(mol_cart, key)._array[0,9:15]
        check_pure = np.dot(block_cart, tfs[2].T)
        error = abs(block_pure - check_pure).max()
        assert error < 1e-5


def gb4_helper(sign0, sign1, sign2, sign3):
    assert abs(sign0) == 1
    assert abs(sign1) == 1
    assert abs(sign2) == 1
    assert abs(sign3) == 1
    max_shell_type = 4
    max_nbasis = get_shell_nbasis(max_shell_type)
    r0 = np.array([ 0.57092,  0.29608, -0.758  ])
    r1 = np.array([ 0.83984,  0.65053,  0.36087])
    r2 = np.array([-0.70841,  0.22864,  0.79589])
    r3 = np.array([-0.62267, -0.83676, -0.75233])
    scales0 = np.ones(15, float)
    scales1 = np.ones(10, float)
    scales2 = np.ones(6, float)
    scales3 = np.ones(6, float)

    gb4oi = GB4ElectronRepulsionIntegralLibInt(max_shell_type)
    assert gb4oi.max_nbasis == max_nbasis
    assert gb4oi.nwork == max_nbasis**4

    gb4oi.reset(sign0*4, sign1*3, sign2*2, sign3*2, r0, r1, r2, r3)
    gb4oi.add(1.0, 5.398, 0.320, 0.123, 0.210, scales0, scales1, scales2, scales3)
    gb4oi.add(0.5, 0.123, 0.210, 1.234, 2.333, scales0, scales1, scales2, scales3)
    gb4oi.add(0.7, 1.234, 2.333, 0.500, 0.500, scales0, scales1, scales2, scales3)
    gb4oi.add(0.3, 0.500, 0.500, 5.398, 0.320, scales0, scales1, scales2, scales3)
    work0 = gb4oi.get_work(15, 10, 6, 6)
    gb4oi.cart_to_pure()
    n0 = get_shell_nbasis(sign0*4)
    n1 = get_shell_nbasis(sign1*3)
    n2 = get_shell_nbasis(sign2*2)
    n3 = get_shell_nbasis(sign3*2)
    work1 = gb4oi.get_work(n0, n1, n2, n3)
    return work0, work1

def test_gb4_electron_repulsion_integral_class_pppp():
    work0, work1 = gb4_helper(1, 1, 1, 1)
    assert work0.shape == work1.shape
    assert abs(work0-work1).max() < 1e-10

def test_gb4_electron_repulsion_integral_class_pppm1():
    work0, work1 = gb4_helper(1, 1, 1, -1)
    work0 = work0.reshape(-1, work0.shape[-1])
    work1 = work1.reshape(-1, work1.shape[-1])
    assert work0.shape[0] == work1.shape[0]
    work0 = np.dot(work0, tfs[2].T)
    assert work0.shape == work1.shape
    assert abs(work0-work1).max() < 1e-10

def test_gb4_electron_repulsion_integral_class_pppm2():
    work0, work1 = gb4_helper(1, 1, 1, -1)
    work0 = np.tensordot(work0, tfs[2], ([3,1]))
    assert work0.shape == work1.shape
    assert abs(work0-work1).max() < 1e-10

def test_gb4_electron_repulsion_integral_class_ppmp():
    work0, work1 = gb4_helper(1, 1, -1, 1)
    work0 = np.tensordot(work0, tfs[2], ([2,1])).transpose(0,1,3,2)
    assert work0.shape == work1.shape
    assert abs(work0-work1).max() < 1e-10

def test_gb4_electron_repulsion_integral_class_pmpp():
    work0, work1 = gb4_helper(1, -1, 1, 1)
    work0 = np.tensordot(work0, tfs[3], ([1,1])).transpose(0,3,1,2)
    assert work0.shape == work1.shape
    assert abs(work0-work1).max() < 1e-10

def test_gb4_electron_repulsion_integral_class_mppp():
    work0, work1 = gb4_helper(-1, 1, 1, 1)
    work0 = np.tensordot(work0, tfs[4], ([0,1])).transpose(3,0,1,2)
    assert work0.shape == work1.shape
    assert abs(work0-work1).max() < 1e-10
