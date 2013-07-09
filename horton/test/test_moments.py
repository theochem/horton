# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


from nose.tools import assert_raises
from horton import *
from horton.test.common import get_pentagon_moments


def test_get_ncart():
    assert get_ncart(0) == 1
    assert get_ncart(1) == 3
    assert get_ncart(2) == 6
    assert get_ncart(3) == 10


def test_get_ncart_cumul():
    assert get_ncart_cumul(0) == 1
    assert get_ncart_cumul(1) == 4
    assert get_ncart_cumul(2) == 10
    assert get_ncart_cumul(3) == 20


def test_get_cartesian_powers():
    lmax = 4
    cartesian_powers = get_cartesian_powers(lmax)
    assert issubclass(cartesian_powers.dtype.type, int)
    assert cartesian_powers.shape == (get_ncart_cumul(lmax), 3)
    assert (cartesian_powers[0] == [0, 0, 0]).all()
    assert (cartesian_powers[1] == [1, 0, 0]).all()
    assert (cartesian_powers[4] == [2, 0, 0]).all()
    assert (cartesian_powers[5] == [1, 1, 0]).all()
    assert (cartesian_powers[6] == [1, 0, 1]).all()
    assert (cartesian_powers[10] == [3, 0, 0]).all()
    assert (cartesian_powers[11] == [2, 1, 0]).all()
    assert (cartesian_powers[19] == [0, 0, 3]).all()
    assert (cartesian_powers[-1] == [0, 0, 4]).all()

    for lmax in xrange(4):
        tmp = get_cartesian_powers(lmax)
        assert tmp.shape == (get_ncart_cumul(lmax), 3)
        assert (tmp == cartesian_powers[:len(tmp)]).all()


def test_rotate_cartesian_moments_mult():
    for mult in 2, 3, 4, 5, 6, 7:
        axis = np.random.normal(0, 1, 3)
        rmat = get_rotation_matrix(axis, 2.0*np.pi/mult)
        m0 = get_pentagon_moments()
        m1 = get_pentagon_moments()
        for i in xrange(mult):
            m1 = rotate_cartesian_moments(m1, rmat)
        assert abs(m0 - m1).max() < 1e-10


def test_rotate_cartesian_moments_general():
    for i in xrange(10):
        rmat = get_rotation_matrix(np.array([1,1,1]), 2.0*np.pi/3)
        m0 = get_pentagon_moments()
        m1 = get_pentagon_moments(rmat)
        m2 = rotate_cartesian_moments(m0, rmat)
        assert abs(m1 - m2).max() < 1e-10


def test_fill_cartesian_polynomials():
    cps = get_cartesian_powers(4)
    output = np.zeros(get_ncart_cumul(4)-1)
    output[:3] = np.random.normal(0, 1, 3)
    assert fill_cartesian_polynomials(output, 0) == -1
    for l in xrange(1, 5):
        assert fill_cartesian_polynomials(output, l) == get_ncart_cumul(l-1)-1
        nrow = get_ncart_cumul(l)-1
        for irow in xrange(nrow):
            px, py, pz = cps[irow+1]
            check = output[0]**px * output[1]**py * output[2]**pz
            assert abs(output[irow] - check) < 1e-10
        assert (output[nrow:] == 0).all()

    with assert_raises(ValueError):
        fill_cartesian_polynomials(output, 5)
    with assert_raises(ValueError):
        fill_cartesian_polynomials(output[:get_ncart_cumul(4)-2], 4)


def test_fill_radial_polynomials():
    lmax = 4
    output = np.zeros(lmax)
    output[0] = 1.13
    assert fill_radial_polynomials(output, lmax) is None
    assert abs(output[1] - output[0]**2) < 1e-10
    assert abs(output[2] - output[0]**3) < 1e-10
    assert abs(output[3] - output[0]**4) < 1e-10
