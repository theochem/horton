# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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


def get_fake_system():
    coordinates = np.array([[0.0, 1.5, 2.3], [-0.1, 1.1, 0.0], [2.0, 1.0, 0.0], [-1.0, 1.0, 1.1]])
    numbers = np.array([1, 1, 2, 2])
    origin = np.array([1.0, 0.0, 1.0])
    grid_rvecs = np.array([[0.15, 0.0, 0.0], [0.0, 0.20, 0.01], [0.01, 0.01, 0.15]])
    shape = np.array([10, 10, 20])
    pbc = np.array([1, 0, 1])
    ugrid = UniformGrid(origin, grid_rvecs, shape, pbc)
    return coordinates, numbers, ugrid


def test_mask_dens():
    coordinates, numbers, ugrid = get_fake_system()
    rho = np.zeros(ugrid.shape)
    scan = np.arange(-2.0, -0.0001, 0.1)
    rho[0,0,:] = 10**scan
    weights = setup_weights(coordinates, numbers, ugrid, dens=(rho, -9, 0.8))
    assert (weights[1,:,:] == 0.0).all()
    assert abs(weights[0,0,:] - np.exp(-0.8*(np.log(rho[0,0,:])-(-9))**2)).max() < 1e-10


def test_mask_near1():
    coordinates, numbers, ugrid = get_fake_system()
    weights = setup_weights(coordinates, numbers, ugrid, near={1: (0.5, 0.5), 2: (1.0, 0.2)})
    assert (weights >= 0.0).all()
    assert (weights <= 1.0).all()
    # find the point close to atom 2 and check that the weight is zero
    grid_cell = ugrid.get_grid_cell()
    i = np.round(grid_cell.to_frac(coordinates[2] - ugrid.origin)).astype(int)
    i[0] = i[0]%10
    i[2] = i[2]%20
    assert weights[i[0], i[1], i[2]] == 0.0


def test_mask_near2():
    coordinates, numbers, ugrid = get_fake_system()
    weights = setup_weights(coordinates, numbers, ugrid, near={1: (0.5, 0.5), 2: (1.0, 0.2)})
    weights1 = setup_weights(coordinates, numbers, ugrid, near={1: (0.5, 0.5)})
    weights2 = setup_weights(coordinates, numbers, ugrid, near={2: (1.0, 0.2)})
    assert abs(weights - weights1*weights2).max() < 1e-10


def test_mask_near3():
    coordinates, numbers, ugrid = get_fake_system()
    weights = setup_weights(coordinates, numbers, ugrid, near={0: (0.5, 0.5)})
    weights1 = setup_weights(coordinates, numbers, ugrid, near={1: (0.5, 0.5)})
    weights2 = setup_weights(coordinates, numbers, ugrid, near={2: (0.5, 0.5)})
    assert abs(weights - weights1*weights2).max() < 1e-10


def test_mask_far():
    coordinates, numbers, ugrid = get_fake_system()
    weights = setup_weights(coordinates, numbers, ugrid, far=(1.0, 0.5))
    assert (weights >= 0.0).all()
    assert (weights <= 1.0).all()
    # find the point close to atom 2 and check that the weight is one
    grid_cell = ugrid.get_grid_cell()
    i = np.round(grid_cell.to_frac(coordinates[2] - ugrid.origin)).astype(int)
    i[0] = i[0]%10
    i[2] = i[2]%20
    assert weights[i[0], i[1], i[2]] == 1.0
