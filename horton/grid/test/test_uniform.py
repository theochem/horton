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

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_uig_integrate_gauss():
    # grid parameters
    spacing = 0.1
    naxis = 81

    # grid setup
    offset = 0.5*spacing*(naxis-1)
    origin = np.zeros(3, float)-offset
    rvecs = np.identity(3)*spacing
    shape = np.array([naxis, naxis, naxis])
    pbc = np.ones(3, int)
    uig = UniformGrid(origin, rvecs, shape, pbc)

    # fill a 3D grid with a gaussian function
    x, y, z = np.indices(shape)
    dsq = (x*spacing-offset)**2 + (y*spacing-offset)**2 + (z*spacing-offset)**2
    data = np.exp(-0.5*dsq)/(2*np.pi)**(1.5)

    # compute the integral and compare with analytical result
    num_result = uig.integrate(data)
    assert abs(num_result - 1.0) < 1e-3


def test_ugrid_writable():
    origin = np.random.normal(0, 1, 3)
    grid_rvecs = np.random.normal(0, 1, (3, 3))
    shape = np.random.randint(10, 20, 3)
    pbc = np.random.randint(0, 2, 3)
    ugrid = UniformGrid(origin, grid_rvecs, shape, pbc)
    # write to attributes and check value
    for i in xrange(10):
        value = np.random.normal(0, 1)
        i, j = np.random.randint(0, 3, 2)
        ugrid.origin[i] = value
        assert ugrid.origin[i] == value
        ugrid.grid_rvecs[i,j] = value
        assert ugrid.grid_rvecs[i,j] == value
        k = np.random.randint(10, 20)
        ugrid.shape[i] = k
        assert ugrid.shape[i] == k
        l = np.random.randint(0, 2)
        ugrid.pbc[j] = l
        assert ugrid.pbc[j] == l


def test_index_wrap():
    assert index_wrap(2, 3) == 2
    assert index_wrap(-5, 10) == 5
    assert index_wrap(5, 10) == 5
    assert index_wrap(15, 10) == 5


def get_simple_test_uig():
    origin = np.array([0.1, -2.0, 3.1])
    rvecs = np.array([
        [1.0, 0.1, 0.2],
        [0.1, 1.1, 0.2],
        [0.0, -0.1, 1.0],
    ])
    #origin = np.array([0.0, 0.0, 0.0])
    #rvecs = np.array([
    #    [1.0, 0.0, 0.0],
    #    [0.0, 1.0, 0.0],
    #    [0.0, 0.0, 1.0],
    #])
    shape = np.array([40, 40, 40])
    pbc = np.array([1, 0, 1])
    return UniformGrid(origin, rvecs, shape, pbc)


def test_get_ranges_rcut1():
    uig = get_simple_test_uig()
    center = np.array([0.1, -2.5, 3.2])
    rb1, re1 = uig.get_ranges_rcut(center, 2.0)
    rb2, re2 = uig.get_grid_cell().get_ranges_rcut(uig.origin - center, 2.0)
    assert rb1[0] == rb2[0]
    assert rb1[1] == 0
    assert rb1[2] == rb2[2]
    assert (re1 == re2).all()


def test_get_ranges_rcut2():
    uig = get_simple_test_uig()
    center = np.array([60.0, 50.0, 60.0])
    rb1, re1 = uig.get_ranges_rcut(center, 2.0)
    rb2, re2 = uig.get_grid_cell().get_ranges_rcut(uig.origin - center, 2.0)
    assert (rb1 == rb2).all()
    assert re1[0] == re2[0]
    assert re1[1] == 40
    assert re1[2] == re2[2]


def test_dist_grid_point():
    uig = get_simple_test_uig()
    assert uig.dist_grid_point(uig.origin, np.array([0, 0, 0])) == 0.0
    assert uig.dist_grid_point(uig.origin, np.array([0, 0, 1])) == (0.1**2+1.0)**0.5
    assert uig.dist_grid_point(uig.origin, np.array([0, 1, 0])) == (0.1**2+1.1**2+0.2**2)**0.5

    center = np.array([0.9, 2.5, 1.6])
    indexes = np.array([6, 3, -2])
    point = uig.origin + uig.get_grid_cell().to_cart(indexes.astype(float))
    assert abs(uig.dist_grid_point(center, np.array([6, 3, -2])) - np.linalg.norm(center - point)) < 1e-10


def test_delta_grid_point():
    uig = get_simple_test_uig()
    assert (uig.delta_grid_point(uig.origin, np.array([0, 0, 0])) == np.array([0.0, 0.0, 0.0])).all()
    assert (uig.delta_grid_point(uig.origin, np.array([0, 0, 1])) == np.array([0.0, -0.1, 1.0])).all()
    assert (uig.delta_grid_point(uig.origin, np.array([0, 1, 0])) == np.array([0.1, 1.1, 0.2])).all()

    center = np.array([0.9, 2.5, 1.6])
    indexes = np.array([6, 3, -2])
    point = uig.origin + uig.get_grid_cell().to_cart(indexes.astype(float))
    assert abs(uig.delta_grid_point(center, np.array([6, 3, -2])) - (point - center)).max() < 1e-10
