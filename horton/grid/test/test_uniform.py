# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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


import numpy as np

from horton import *
from horton.grid.test.common import *
from horton.test.common import get_random_cell


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
    uig = UniformIntGrid(origin, rvecs, shape, pbc)

    # fill a 3D grid with a gaussian function
    x, y, z = np.indices(shape)
    dsq = (x*spacing-offset)**2 + (y*spacing-offset)**2 + (z*spacing-offset)**2
    data = np.exp(-0.5*dsq)/(2*np.pi)**(1.5)

    # compute the integral and compare with analytical result
    num_result = uig.integrate(data)
    assert abs(num_result - 1.0) < 1e-3


def test_index_wrap():
    assert index_wrap(2, 3) == 2
    assert index_wrap(-5, 10) == 5
    assert index_wrap(5, 10) == 5
    assert index_wrap(15, 10) == 5


def test_uig_eval_spline_simple1():
    cs = get_cosine_spline()
    offset = 5.0
    spacing = 0.5
    uig = UniformIntGrid(np.array([-offset, -offset, -offset]), np.identity(3)*spacing, np.array([21,21,21]), np.array([1, 1, 1]))

    x, y, z = np.indices(uig.shape)
    d = np.sqrt((x*spacing-offset)**2 + (y*spacing-offset)**2 + (z*spacing-offset)**2).ravel()
    data1 = cs(d)
    data2 = np.zeros(uig.shape)
    uig.eval_spline(cs, np.zeros(3), data2)
    data2 = data2.ravel()

    assert abs(data1 - data2).max() == 0.0


def test_uig_eval_spline_simple2():
    cs = get_cosine_spline()
    offset = 1.0
    spacing = 0.1
    uig = UniformIntGrid(np.array([-offset, -offset, -offset]), np.identity(3)*spacing, np.array([21,21,21]), np.array([1, 1, 1]))

    x, y, z = np.indices(uig.shape)
    data1 = 0
    for ix in xrange(-1,2):
        for iy in xrange(-1,2):
            for iz in xrange(-1,2):
                d = np.sqrt((x*spacing-offset+ix*2.1)**2 + (y*spacing-offset+iy*2.1)**2 + (z*spacing-offset+iz*2.1)**2).ravel()
                data1 += cs(d)

    data2 = np.zeros(uig.shape)
    uig.eval_spline(cs, np.zeros(3), data2)
    data2 = data2.ravel()

    assert abs(data1 - data2).max() < 1e-12


def test_uig_eval_spline_with_integration():
    cs = get_cosine_spline()

    # Different integration grids
    pbc = np.ones(3, int)
    uigs = [
        UniformIntGrid(np.array([-1.0, -1.0, -1.0]), np.identity(3)*0.1, np.array([21,21,21]), pbc),
        UniformIntGrid(np.array([-1.0, -1.0, -1.0]), np.array([[0.1, 0.1, 0.0], [0.1, 0.0, 0.1], [0.0, 0.1, 0.1]]), np.array([21,21,21]), pbc),
        UniformIntGrid(np.array([0.0, 0.0, 0.0]), np.identity(3)*0.1, np.array([21,21,21]), pbc),
        UniformIntGrid(np.array([0.0, 0.0, 0.0]), np.array([[0.1, 0.1, 0.0], [0.1, 0.0, 0.1], [0.0, 0.1, 0.1]]), np.array([21,21,21]), pbc),
        UniformIntGrid(np.array([0.0, 0.0, 0.0]), np.array([[0.1471, 0.0745, -0.010], [0.0315, -0.1403, 0.1], [0.0014, 0.147, 0.0]]), np.array([21,21,21]), pbc),
    ]

    for uig in uigs:
        # fill a 3D grid with the cosine test function
        data = np.zeros(uig.shape)
        uig.eval_spline(cs, np.zeros(3), data)

        # test the integral
        expected = 4*np.pi**2*(np.pi**2/3-2)
        assert abs(uig.integrate(data) - expected) < 1e-3


def test_uig_eval_spline_3d_random():
    cs = get_cosine_spline()

    for i in xrange(10):
        origin = np.random.uniform(-1, 1, 3)
        grid_cell = get_random_cell(0.3, 3)
        shape = np.random.randint(10, 20, 3)
        pbc = np.array([1, 1, 1])
        uig = UniformIntGrid(origin, grid_cell.rvecs, shape, pbc)

        rvecs = grid_cell.rvecs*uig.shape.reshape(-1,1)

        output1 = np.zeros(uig.shape)
        center1 = np.random.uniform(-3, 3, 3)
        uig.eval_spline(cs, center1, output1)
        output2 = np.zeros(uig.shape)
        center2 = center1 + np.dot(np.random.randint(-3, 3, 3), rvecs)
        uig.eval_spline(cs, center2, output2)

        assert abs(output1 - output2).max() < 1e-10


def test_uig_eval_spline_2d_random():
    cs = get_cosine_spline()

    for i in xrange(15):
        origin = np.random.uniform(-1, 1, 3)
        grid_cell = get_random_cell(0.2, 3)
        shape = np.random.randint(10, 20, 3)
        pbc = np.array([1, 1, 1])
        pbc[np.random.randint(3)] = 0
        uig = UniformIntGrid(origin, grid_cell.rvecs, shape, pbc)

        tmp = uig.shape*pbc
        rvecs = grid_cell.rvecs*tmp.reshape(-1,1)

        output1 = np.zeros(uig.shape)
        center1 = np.random.uniform(0, 1, 3)
        uig.eval_spline(cs, center1, output1)
        output2 = np.zeros(uig.shape)
        center2 = center1 + np.dot(np.random.randint(-3, 3, 3), rvecs)
        uig.eval_spline(cs, center2, output2)

        assert abs(output1 - output2).max() < 1e-10


def test_uig_eval_spline_0d_random():
    cs = get_cosine_spline()

    for i in xrange(15):
        origin = np.random.uniform(-1, 1, 3)
        grid_cell = get_random_cell(0.2, 3)
        shape = np.random.randint(10, 20, 3)
        pbc = np.array([0, 0, 0])
        uig = UniformIntGrid(origin, grid_cell.rvecs, shape, pbc)

        center = np.random.uniform(0, 2, 3)
        output1 = np.zeros(uig.shape)
        uig.eval_spline(cs, center, output1)

        x, y, z = np.indices(shape)
        x = x.ravel()
        y = y.ravel()
        z = z.ravel()
        rvecs = grid_cell.rvecs
        points = np.outer(x, rvecs[0]) + np.outer(y, rvecs[1]) + np.outer(z, rvecs[2])
        points += origin
        distances = np.zeros(len(points))
        grid_distances(points, center, distances)
        output2 = cs(distances)
        output2.shape = shape

        assert abs(output1 - output2).max() < 1e-10


def test_uig_eval_spline_1d_random():
    cs = get_cosine_spline()

    for i in xrange(15):
        origin = np.random.uniform(-1, 1, 3)
        grid_cell = get_random_cell(0.2, 3)
        shape = np.random.randint(10, 20, 3)
        pbc = np.array([0, 0, 0])
        pbc[np.random.randint(3)] = 1
        uig = UniformIntGrid(origin, grid_cell.rvecs, shape, pbc)

        tmp = uig.shape*pbc
        rvecs = grid_cell.rvecs*tmp.reshape(-1,1)

        output1 = np.zeros(uig.shape)
        center1 = np.random.uniform(0, 1, 3)
        uig.eval_spline(cs, center1, output1)
        output2 = np.zeros(uig.shape)
        center2 = center1 + np.dot(np.random.randint(-3, 3, 3), rvecs)
        uig.eval_spline(cs, center2, output2)

        assert abs(output1 - output2).max() < 1e-10


def test_uig_eval_spline_add_random():
    cs = get_cosine_spline()

    for i in xrange(20):
        origin = np.random.uniform(-1, 1, 3)
        grid_cell = get_random_cell(0.2, 3)
        shape = np.random.randint(10, 20, 3)
        pbc = np.random.randint(0, 2, 3).astype(int)
        uig = UniformIntGrid(origin, grid_cell.rvecs, shape, pbc)

        output1 = np.zeros(uig.shape)
        center1 = np.random.uniform(0, 1, 3)
        uig.eval_spline(cs, center1, output1)

        output2 = np.zeros(uig.shape)
        center2 = np.random.uniform(0, 1, 3)
        uig.eval_spline(cs, center2, output2)

        output3 = np.zeros(uig.shape)
        uig.eval_spline(cs, center1, output3)
        uig.eval_spline(cs, center2, output3)

        assert abs(output1 + output2 - output3).max() < 1e-10


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
    pbc = np.ones(3, int)
    return UniformIntGrid(origin, rvecs, shape, pbc)


def test_weight_corrections():
    from horton.cpart.test.common import get_fake_co

    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    funcs = [
        (sys.coordinates[0], [
            #proatomdb.get_spline(6, +1),
            proatomdb.get_spline(6,  0),
            #proatomdb.get_spline(6, -1),
        ]),
        (sys.coordinates[1], [
            #proatomdb.get_spline(8, +1),
            proatomdb.get_spline(8,  0),
            #proatomdb.get_spline(8, -1),
        ]),
    ]
    weights = ui_grid.compute_weight_corrections(funcs)
    assert (weights != 1.0).any()
    assert (weights == 1.0).any()

    assert abs(ui_grid.integrate(mol_dens, weights)-14.0) < 6e-3
    assert abs(ui_grid.integrate(mol_dens)-14.0) > 5e-2
