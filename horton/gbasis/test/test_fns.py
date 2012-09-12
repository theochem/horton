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


import numpy as np

from horton import *


def test_exceptions():
    try:
        grid_fn = GB1GridFn(-1)
        assert False
    except ValueError:
        pass

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])

    try:
        grid_fn = GB1GridFn(2)
        grid_fn.reset(-3, center, point)
        assert False
    except ValueError:
        pass

    try:
        grid_fn = GB1GridFn(2)
        grid_fn.reset(3, center, point)
        assert False
    except ValueError:
        pass


def test_grid_fn_s():
    grid_fn = GB1GridFn(0)
    assert grid_fn.nwork == 1
    assert grid_fn.max_shell_type == 0
    assert grid_fn.max_nbasis == 1

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])
    grid_fn.reset(0, center, point)
    assert grid_fn.shell_type0 == 0

    coeff = 0.3
    alpha = 0.5
    scale0 = 0.7
    grid_fn.add(coeff, alpha, np.array([scale0]))
    work = grid_fn.get_work(1)
    assert work.shape == (1,)

    dsq = np.linalg.norm(center - point)**2
    assert abs(work[0] -  scale0*coeff*np.exp(-alpha*dsq)) < 1e-10

def test_grid_fn_p():
    grid_fn = GB1GridFn(1)
    assert grid_fn.nwork == 3
    assert grid_fn.max_shell_type == 1
    assert grid_fn.max_nbasis == 3

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])
    grid_fn.reset(1, center, point)
    assert grid_fn.shell_type0 == 1

    coeff = 0.3
    alpha = 0.5
    scales0 = np.array([0.1, 0.2, 0.7])
    grid_fn.add(coeff, alpha, scales0)
    work = grid_fn.get_work(3)
    assert work.shape == (3,)

    d = point - center
    dsq = np.linalg.norm(d)**2
    for i in xrange(3):
        assert abs(work[i] -  scales0[i]*coeff*np.exp(-alpha*dsq)*d[i]) < 1e-10


def test_grid_fn_p_contraction():
    grid_fn = GB1GridFn(1)
    assert grid_fn.nwork == 3
    assert grid_fn.max_shell_type == 1
    assert grid_fn.max_nbasis == 3

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])
    grid_fn.reset(1, center, point)
    assert grid_fn.shell_type0 == 1

    scales0 = np.array([0.1, 0.2, 0.7])
    coeff0 = 0.3
    alpha0 = 0.5
    grid_fn.add(coeff0, alpha0, scales0)
    coeff1 = 0.7
    alpha1 = 0.02
    grid_fn.add(coeff1, alpha1, scales0)

    work = grid_fn.get_work(3)
    assert work.shape == (3,)

    d = point - center
    dsq = np.linalg.norm(d)**2
    for i in xrange(3):
        expected = scales0[i]*coeff0*np.exp(-alpha0*dsq)*d[i] + \
                   scales0[i]*coeff1*np.exp(-alpha1*dsq)*d[i]
        assert abs(work[i] - expected) < 1e-10


def test_grid_fn_d_contraction():
    grid_fn = GB1GridFn(3)
    assert grid_fn.nwork == 10
    assert grid_fn.max_shell_type == 3
    assert grid_fn.max_nbasis == 10

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])
    grid_fn.reset(-2, center, point)
    assert grid_fn.shell_type0 == -2

    scales0 = np.array([0.1, 0.2, 0.7, 0.6, 0.3, 0.8])
    coeff0 = 0.3
    alpha0 = 0.5
    grid_fn.add(coeff0, alpha0, scales0)
    coeff1 = 0.7
    alpha1 = 0.02
    grid_fn.add(coeff1, alpha1, scales0)

    work_cart = grid_fn.get_work(6)
    assert work_cart.shape == (6,)

    d = point - center
    dsq = np.linalg.norm(d)**2

    expected = scales0[0]*coeff0*np.exp(-alpha0*dsq)*d[0]*d[0] + \
               scales0[0]*coeff1*np.exp(-alpha1*dsq)*d[0]*d[0]
    assert abs(work_cart[0] - expected) < 1e-10 # xx

    expected = scales0[3]*coeff0*np.exp(-alpha0*dsq)*d[1]*d[1] + \
               scales0[3]*coeff1*np.exp(-alpha1*dsq)*d[1]*d[1]
    assert abs(work_cart[3] - expected) < 1e-10 # yy

    expected = scales0[4]*coeff0*np.exp(-alpha0*dsq)*d[1]*d[2] + \
               scales0[4]*coeff1*np.exp(-alpha1*dsq)*d[1]*d[2]
    assert abs(work_cart[4] - expected) < 1e-10 # yz

    grid_fn.cart_to_pure()
    work_pure = grid_fn.get_work(5)
    assert work_pure.shape == (5,)

    from test_cartpure import tfs
    assert abs(work_pure - np.dot(tfs[2], work_cart)).max() < 1e-10
