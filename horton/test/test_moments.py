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

from horton.test.common import get_pentagon_moments, get_point_moments


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


def test_rotate_cartesian_moments_pentagon_mult():
    for mult in 2, 3, 4, 5:
        axis = np.random.normal(0, 1, 3)
        rmat = get_rotation_matrix(axis, np.pi/mult)
        m0 = get_pentagon_moments()
        m1 = get_pentagon_moments()
        for i in xrange(mult):
            m1 = rotate_cartesian_moments_all(rmat, m1)
            m1 = rotate_cartesian_moments_all(rmat, m1)
        assert abs(m0 - m1).max() < 1e-10


def test_rotate_cartesian_moments_pentagon_general():
    for i in xrange(10):
        rmat = get_random_rotation()
        m0 = get_pentagon_moments()
        m1 = get_pentagon_moments(rmat)
        m2 = rotate_cartesian_moments_all(rmat, m0)
        assert abs(m1 - m2).max() < 1e-10


def test_rotate_cartesian_moments_random_mult():
    for i in xrange(10):
        coordinates = np.random.normal(0, 1, (10, 3))
        for mult in 2, 3, 4, 5:
            axis = np.random.normal(0, 1, 3)
            rmat = get_rotation_matrix(axis, np.pi/mult)
            m0 = get_point_moments(coordinates)
            m1 = get_point_moments(coordinates)
            for i in xrange(mult):
                m1 = rotate_cartesian_moments_all(rmat, m1)
                m1 = rotate_cartesian_moments_all(rmat, m1)
            assert abs(m0 - m1).max() < 1e-10


def test_rotate_cartesian_moments_random_general():
    for i in xrange(10):
        coordinates = np.random.normal(0, 1, (10, 3))
        rmat = get_random_rotation()
        m0 = get_point_moments(coordinates)
        m1 = get_point_moments(coordinates, rmat)
        m2 = rotate_cartesian_moments_all(rmat, m0)
        assert abs(m1 - m2).max() < 1e-10


def test_rotate_simple_moments():
    rmat = get_rotation_matrix(np.array([0, 0, 1]), np.pi/4)
    m0 = np.array([1, 0, 0])
    m1 = rotate_cartesian_multipole(rmat, m0, 'moments')
    assert abs(m1 - [0.5**0.5, 0.5**0.5, 0.0]).max() < 1e-10
    m0 = np.array([1, 0, 0, 0, 0, 0])
    m1 = rotate_cartesian_multipole(rmat, m0, 'moments')
    assert abs(m1 - [0.5, 0.5, 0, 0.5, 0.0, 0.0]).max() < 1e-10


def test_rotate_simple_coeffs():
    rmat = get_rotation_matrix(np.array([0, 0, 1]), np.pi/4)
    m0 = np.array([1, 0, 0])
    m1 = rotate_cartesian_multipole(rmat, m0, 'coeffs')
    assert abs(m1 - [0.5**0.5, 0.5**0.5, 0.0]).max() < 1e-10
    m0 = np.array([1, 0, 0, 0, 0, 0])
    m1 = rotate_cartesian_multipole(rmat, m0, 'coeffs')
    assert abs(m1 - [0.5, 1.0, 0, 0.5, 0.0, 0.0]).max() < 1e-10


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


def test_get_npure():
    assert get_npure(0) == 1
    assert get_npure(1) == 3
    assert get_npure(2) == 5
    assert get_npure(3) == 7


def test_get_npure_cumul():
    assert get_npure_cumul(0) == 1
    assert get_npure_cumul(1) == 4
    assert get_npure_cumul(2) == 9
    assert get_npure_cumul(3) == 16


def test_fill_pure_polynomials():
    output = np.zeros(get_npure_cumul(4)-1)
    #output[:3] = np.random.normal(0, 1, 3)
    output[:3] = [0.35, 0.12, 0.55]
    z = output[0]
    x = output[1]
    y = output[2]
    r = np.linalg.norm(output[:3])
    expected = np.array([
        # l = 1
        z, x, y,
        # l = 2
        (3.0*z**2 - r**2)/2.0,
        3.0**0.5*x*z,
        3.0**0.5*y*z,
        3.0**0.5*(x*x-y*y)/2.0,
        3.0**0.5*x*y,
        # l = 3
        (5.0*z**3 - 3.0*r**2.0*z)/2,
        x*(15.0*z**2 - 3.0*r**2)/2.0/6.0**0.5,
        y*(15.0*z**2 - 3.0*r**2)/2.0/6.0**0.5,
        z*(x**2-y**2)*15.0**0.5/2,
        x*y*z*15.0**0.5,
        10.0**0.5/4*x*(x**2-3.0*y**2),
        10.0**0.5/4*y*(3.0*x**2-y**2),
        # l = 4
        (3.0*r**4 - 30.0*r**2*z**2 + 35.0*z**4)/8,
        x*(35*z**3 - 15.0*r**2*z)/2.0/10.0**0.5,
        y*(35*z**3 - 15.0*r**2*z)/2.0/10.0**0.5,
        (105.0*z**2-15.0*r**2)*(x**2-y**2)/2.0/30.0*5.0**0.5,
        x*y*(105.0*z**2-15.0*r**2)/2.0/15.0*5**0.5,
        z*x*(x**2 - 3*y**2)*70**0.5/4.0,
        z*y*(3*x**2 - y**2)*70**0.5/4.0,
        (x**4 - 6.0*x**2*y**2 + y**4)*35**0.5/8,
        x*y*(x**2-y**2)*35**0.5/2,
    ])
    assert fill_pure_polynomials(output, 0) == -1
    for l in xrange(1, 5):
        assert fill_pure_polynomials(output, l) == get_npure_cumul(l-1)-1
        nrow = get_npure_cumul(l)-1
        #print l, nrow
        for irow in xrange(nrow):
            #print irow, output[irow], expected[irow]
            assert abs(output[irow] - expected[irow]) < 1e-10
        #print
        assert (output[nrow:] == 0).all()

    with assert_raises(ValueError):
        fill_cartesian_polynomials(output, 5)
    with assert_raises(ValueError):
        fill_cartesian_polynomials(output[:get_ncart_cumul(4)-2], 4)


def test_fill_pure_polynomials_array():
    lmax = 4
    npoint = 10
    work = np.zeros( (npoint, (lmax+1)**2-1) )
    work[:,:3] = np.random.normal(0, 1, (npoint, 3))
    fill_pure_polynomials(work, lmax)
    for row in work:
        tmp = row.copy()
        tmp[3:] = 0.0
        fill_pure_polynomials(tmp, lmax)
        assert abs(row[3:] - tmp[3:]).max() < 1e-15


def test_ortho_and_norm_pure():
    for radius in 0.5, 1.0, 1.8797:
        nll = 170
        points = np.zeros((nll, 3))
        weights = np.zeros(nll)
        lebedev_laikov_sphere(points, weights)
        points *= radius
        weights *= 4*np.pi

        lmax = 5
        work = np.zeros((nll, (lmax+1)**2-1))
        work[:,0] = points[:,2]
        work[:,1] = points[:,0]
        work[:,2] = points[:,1]
        fill_pure_polynomials(work, lmax)

        # see if we reproduce racah's norm with the pure polynomials
        counter = 0
        for l in xrange(1, lmax+1):
            for m in xrange(-l, l+1):
                num_norm = np.dot(weights, work[:,counter]*work[:,counter])
                norm = 4*np.pi/(2*l+1)*radius**(2*l)
                assert abs(num_norm - norm) < 1e-10
                counter += 1

        # also test orthogonality
        for counter0 in xrange(counter):
            for counter1 in xrange(counter0):
                tmp = np.dot(weights, work[:,counter0]*work[:,counter1])
                assert abs(tmp) < 1e-10


def test_fill_radial_polynomials():
    lmax = 4
    output = np.zeros(lmax)
    output[0] = 1.13
    assert fill_radial_polynomials(output, lmax) is None
    assert abs(output[1] - output[0]**2) < 1e-10
    assert abs(output[2] - output[0]**3) < 1e-10
    assert abs(output[3] - output[0]**4) < 1e-10
