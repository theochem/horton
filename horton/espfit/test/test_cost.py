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
from horton.test.common import get_random_cell


def get_random_esp_cost_cube3d_args():
    coordinates = np.random.normal(0, 1, (5, 3))
    numbers = np.ones(5, int)
    origin = np.random.uniform(-3, 3, 3)
    rvecs = np.diag(np.random.uniform(2.0, 3.0, 3))
    rvecs += np.random.uniform(-0.1, 0.1, (3, 3))
    grid_cell = Cell(rvecs)
    shape = np.array([5, 5, 5])
    pbc_active = np.array([1, 1, 1])
    vref = np.random.normal(0, 1, shape)
    weights = np.random.uniform(0, 1, shape)
    return coordinates, numbers, origin, grid_cell, shape, pbc_active, vref, weights


def check_costs(costs):
    assert abs(costs[0]._A).max() > 1e-3
    assert abs(costs[0]._B).max() > 1e-3
    for i in xrange(9):
        assert abs(costs[i]._A - costs[i+1]._A).max() < 1e-10
        assert abs(costs[i]._B - costs[i+1]._B).max() < 1e-10
        assert abs(costs[i]._C - costs[i+1]._C) < 1e-10


def test_esp_cost_cube3d_invariance_origin():
    # Some parameters
    coordinates, numbers, origin, grid_cell, shape, pbc_active, vref, weights = \
        get_random_esp_cost_cube3d_args()
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        shift = np.random.uniform(-3, 3, 3)
        sys = System(coordinates+shift, np.ones(5, int))
        grid = UniformIntGrid(shift, grid_cell, shape, pbc_active)
        cost = ESPCost(sys, grid, vref, weights)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)


def test_esp_cost_cube3d_invariance_rotation():
    # Some parameters
    coordinates, numbers, origin, grid_cell, shape, pbc_active, vref, weights = \
        get_random_esp_cost_cube3d_args()
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        A = np.random.normal(0, 1, (3, 3))
        A = 0.5*(A+A.T)
        evals, evecs = np.linalg.eigh(A)
        del A
        del evals
        rvecs = np.dot(grid_cell.rvecs, evecs)
        new_grid_cell = Cell(rvecs)
        new_coordinates = np.dot(coordinates-origin, evecs)+origin

        sys = System(new_coordinates, np.ones(5, int))
        grid = UniformIntGrid(origin, new_grid_cell, shape, pbc_active)
        cost = ESPCost(sys, grid, vref, weights)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)


def test_esp_cost_cube3d_invariance_images():
    # Some parameters
    coordinates, numbers, origin, grid_cell, shape, pbc_active, vref, weights = \
        get_random_esp_cost_cube3d_args()
    grid = UniformIntGrid(origin, grid_cell, shape, pbc_active)
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        rvecs = grid.get_cell().rvecs
        new_coordinates = coordinates.copy()
        for j in xrange(len(coordinates)):
            new_coordinates[j] += np.dot(np.random.randint(-3, 4, 3), rvecs)

        sys = System(new_coordinates, np.ones(5, int))
        cost = ESPCost(sys, grid, vref, weights)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)


def test_esp_cost_cube3d_invariance_rcut():
    # Some parameters
    coordinates, numbers, origin, grid_cell, shape, pbc_active, vref, weights = \
        get_random_esp_cost_cube3d_args()
    grid = UniformIntGrid(origin, grid_cell, shape, pbc_active)
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        sys = System(coordinates, np.ones(5, int))
        cost = ESPCost(sys, grid, vref, weights, rcut=np.random.uniform(10, 30), alpha_scale=4.5, gcut_scale=1.5)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)
