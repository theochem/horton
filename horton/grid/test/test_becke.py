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


def test_becke_sum2_one():
    npoint = 100
    points = np.random.uniform(-5, 5, (npoint, 3))
    weights = np.random.uniform(0, 1, npoint)
    grid = IntGrid(points, weights)

    weights0 = np.ones(npoint, float)
    weights1 = np.ones(npoint, float)

    radii = np.array([0.5, 0.8])
    centers = np.array([[1.2, 2.3, 0.1], [-0.4, 0.0, -2.2]])
    becke_helper_atom(points, weights0, radii, centers, 0, 3)
    becke_helper_atom(points, weights1, radii, centers, 1, 3)

    assert abs(weights0+weights1 - 1).max() < 1e-10


def test_becke_sum3_one():
    npoint = 100
    points = np.random.uniform(-5, 5, (npoint, 3))
    weights = np.random.uniform(0, 1, npoint)
    grid = IntGrid(points, weights)

    weights0 = np.ones(npoint, float)
    weights1 = np.ones(npoint, float)
    weights2 = np.ones(npoint, float)

    radii = np.array([0.5, 0.8, 5.0])
    centers = np.array([[1.2, 2.3, 0.1], [-0.4, 0.0, -2.2], [2.2, -1.5, 0.0]])
    becke_helper_atom(points, weights0, radii, centers, 0, 3)
    becke_helper_atom(points, weights1, radii, centers, 1, 3)
    becke_helper_atom(points, weights2, radii, centers, 2, 3)

    assert abs(weights0+weights1+weights2 - 1).max() < 1e-10


def test_becke_special_points():
    radii = np.array([0.5, 0.8, 5.0])
    centers = np.array([[1.2, 2.3, 0.1], [-0.4, 0.0, -2.2], [2.2, -1.5, 0.0]])

    weights = np.ones(3, float)
    becke_helper_atom(centers, weights, radii, centers, 0, 3)
    assert abs(weights[0] - 1.0) < 1e-10
    assert abs(weights[1]) < 1e-10
    assert abs(weights[2]) < 1e-10

    weights = np.ones(3, float)
    becke_helper_atom(centers, weights, radii, centers, 1, 3)
    assert abs(weights[0]) < 1e-10
    assert abs(weights[1] - 1.0) < 1e-10
    assert abs(weights[2]) < 1e-10

    weights = np.ones(3, float)
    becke_helper_atom(centers, weights, radii, centers, 2, 3)
    assert abs(weights[0]) < 1e-10
    assert abs(weights[1]) < 1e-10
    assert abs(weights[2] - 1.0) < 1e-10
