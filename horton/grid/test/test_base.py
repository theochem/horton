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


import numpy as np
from horton import *
from horton.grid.test.common import get_cosine_spline
from horton.test.common import get_random_cell


def test_grid_integrate():
    npoint = 10
    grid = IntGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    pot = np.random.normal(0, 1, npoint)
    dens = np.random.normal(0, 1, npoint)

    # three
    int1 = grid.integrate(pot, dens)
    int2 = (grid.weights*pot*dens).sum()
    assert abs(int1 - int2) < 1e-10

    # two
    int1 = grid.integrate(pot)
    int2 = (grid.weights*pot).sum()
    assert abs(int1 - int2) < 1e-10

    # one
    int1 = grid.integrate()
    int2 = (grid.weights).sum()
    assert abs(int1 - int2) < 1e-10


def test_grid_integrate_moments():
    npoint = 10
    grid = IntGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    dens = np.random.normal(0, 1, npoint)

    center = np.random.normal(0, 1, 3)

    int1 = grid.integrate(dens, center=center, nx=1)
    x = grid.points[:,0]-center[0]
    int2 = (grid.weights*dens*x).sum()
    assert abs(int1 - int2) < 1e-10

    int1 = grid.integrate(dens, center=center, ny=1)
    y = grid.points[:,1]-center[1]
    int2 = (grid.weights*dens*y).sum()
    assert abs(int1 - int2) < 1e-10

    int1 = grid.integrate(dens, center=center, nz=1)
    z = grid.points[:,2]-center[2]
    int2 = (grid.weights*dens*z).sum()
    assert abs(int1 - int2) < 1e-10

    int1 = grid.integrate(dens, center=center, nr=1)
    r = np.sqrt(sum([(grid.points[:,i]-center[i])**2 for i in xrange(3)]))
    int2 = (grid.weights*dens*r).sum()
    assert abs(int1 - int2) < 1e-10

    int1 = grid.integrate(dens, center=center, nx=2, ny=3, nz=4, nr=1)
    int2 = (grid.weights*dens*x*x*y*y*y*z*z*z*z*r).sum()
    assert abs(int1 - int2) < 1e-10


def test_dot_multi():
    npoint = 10
    pot = np.random.normal(0, 1, npoint)
    dens = np.random.normal(0, 1, npoint)

    # two
    dot1 = np.dot(pot, dens)
    dot2 = dot_multi(pot, dens)
    assert abs(dot1 - dot2) < 1e-10

    # one
    dot1 = pot.sum()
    dot2 = dot_multi(pot)
    assert abs(dot1 - dot2) < 1e-10

    # zero
    dot2 = dot_multi()
    assert dot2 == 0.0


def test_eval_spline_grid_simplest():
    npoint = 10
    boxsize = 2.0
    points = np.random.normal(0, boxsize, (npoint,3))
    #points = np.hstack([np.arange(0.0, 0.9001, 0.1).reshape(-1,1)]*3)
    g = IntGrid(points, np.random.normal(0, 1.0, npoint))
    cs = get_cosine_spline()
    cell = Cell(np.identity(3, float)*boxsize)

    output1 = np.zeros(npoint)
    g.eval_spline(cs, np.array([0.0, 0.0, 0.0]), output1, cell)
    output2 = np.zeros(npoint)
    g.eval_spline(cs, np.array([2.0, 0.0, 0.0]), output2, cell)
    assert abs(output1 - output2).max() < 1e-10


def test_eval_spline_grid_3d_random():
    npoint = 10
    for i in xrange(10):
        cell = get_random_cell(1.0, 3)
        rvecs = cell.rvecs
        points = np.dot(np.random.normal(-2, 3, (npoint,3)), rvecs)
        g = IntGrid(points, np.random.normal(0, 1.0, npoint))
        cs = get_cosine_spline()

        output1 = np.zeros(npoint)
        center1 = np.random.uniform(-10, 10, 3)
        g.eval_spline(cs, center1, output1, cell)
        output2 = np.zeros(npoint)
        center2 = center1 + np.dot(np.random.randint(-3, 3, 3), rvecs)
        g.eval_spline(cs, center2, output2, cell)

        assert abs(output1 - output2).max() < 1e-10


def test_eval_spline_grid_2d_random():
    npoint = 10
    cs = get_cosine_spline()

    for i in xrange(10):
        cell = get_random_cell(1.0, 2)
        rvecs = cell.rvecs
        points = np.dot(np.random.normal(-2, 3, (npoint,2)), rvecs) + np.random.normal(-3, 3, (npoint,3))
        g = IntGrid(points, np.random.normal(0, 1.0, npoint))

        output1 = np.zeros(npoint)
        center1 = np.random.uniform(-3, 3, 3)
        g.eval_spline(cs, center1, output1, cell)
        output2 = np.zeros(npoint)
        center2 = center1 + np.dot(np.random.randint(-3, 3, 2), rvecs)
        g.eval_spline(cs, center2, output2, cell)

        assert abs(output1 - output2).max() < 1e-10


def test_eval_spline_grid_1d_random():
    npoint = 10
    cs = get_cosine_spline()

    for i in xrange(10):
        cell = get_random_cell(1.0, 1)
        rvecs = cell.rvecs
        points = np.random.normal(-3, 3, (npoint,3))
        g = IntGrid(points, np.random.normal(0, 1.0, npoint))

        output1 = np.zeros(npoint)
        center1 = np.random.uniform(-3, 3, 3)
        g.eval_spline(cs, center1, output1, cell)
        output2 = np.zeros(npoint)
        center2 = center1 + np.dot(np.random.randint(-3, 3, 1), rvecs)
        g.eval_spline(cs, center2, output2, cell)

        assert abs(output1 - output2).max() < 1e-10


def test_eval_spline_grid_0d_random():
    npoint = 10
    cs = get_cosine_spline()

    for i in xrange(10):
        cell = Cell(None)
        points = np.random.normal(-1, 1, (npoint,3))
        g = IntGrid(points, np.random.normal(0, 1.0, npoint))

        center = np.random.uniform(-1, 1, 3)

        output1 = np.zeros(npoint)
        g.eval_spline(cs, center, output1, cell)

        distances = np.zeros(npoint)
        g.distances(center, distances)
        output2 = cs(distances)

        assert abs(output1 - output2).max() < 1e-10


def test_eval_spline_grid_add_random():
    npoint = 10
    cs = get_cosine_spline()

    for i in xrange(10):
        cell = get_random_cell(1.0, np.random.randint(4))
        points = np.random.normal(-2, 3, (npoint,3))
        g = IntGrid(points, np.random.normal(0, 1.0, npoint))

        output1 = np.zeros(npoint)
        center1 = np.random.uniform(-2, 2, 3)
        g.eval_spline(cs, center1, output1, cell)

        output2 = np.zeros(npoint)
        center2 = np.random.uniform(-2, 2, 3)
        g.eval_spline(cs, center2, output2, cell)

        output3 = np.zeros(npoint)
        g.eval_spline(cs, center1, output3, cell)
        g.eval_spline(cs, center2, output3, cell)

        assert abs(output1 + output2 - output3).max() < 1e-10
