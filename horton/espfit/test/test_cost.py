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
#pylint: skip-file


import numpy as np, h5py as h5
from horton import *
from horton.test.common import check_delta


def get_random_esp_cost_cube3d_args():
    coordinates = np.random.normal(0, 1, (5, 3))
    numbers = np.ones(5, int)
    origin = np.random.uniform(-3, 3, 3)
    grid_rvecs = np.diag(np.random.uniform(2.0, 3.0, 3))
    grid_rvecs += np.random.uniform(-0.1, 0.1, (3, 3))
    shape = np.array([5, 5, 5])
    pbc = np.array([1, 1, 1])
    vref = np.random.normal(0, 1, shape)
    weights = np.random.uniform(0, 1, shape)
    return coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights


def get_random_esp_cost_cube0d_args():
    coordinates = np.random.normal(0, 1, (5, 3))
    numbers = np.ones(5, int)
    origin = np.random.uniform(-3, 3, 3)
    grid_rvecs = np.diag(np.random.uniform(2.0, 3.0, 3))
    grid_rvecs += np.random.uniform(-0.1, 0.1, (3, 3))
    shape = np.array([5, 5, 5])
    pbc = np.array([0, 0, 0])
    vref = np.random.normal(0, 1, shape)
    weights = np.random.uniform(0, 1, shape)
    return coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights


def get_random_esp_cost_cube3d():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
       get_random_esp_cost_cube3d_args()
    grid = UniformGrid(origin, grid_rvecs, shape, pbc)
    return ESPCost.from_grid_data(coordinates, grid, vref, weights)


def get_random_esp_cost_cube0d():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
       get_random_esp_cost_cube0d_args()
    grid = UniformGrid(origin, grid_rvecs, shape, pbc)
    return ESPCost.from_grid_data(coordinates, grid, vref, weights)


def check_costs(costs, eps0=1e-3, eps1=1e-9):
    assert abs(costs[0]._A).max() > eps0
    assert abs(costs[0]._B).max() > eps0
    for i in xrange(9):
        assert abs(costs[i]._A - costs[i+1]._A).max() < eps1
        assert abs(costs[i]._B - costs[i+1]._B).max() < eps1
        assert abs(costs[i]._C - costs[i+1]._C) < eps1


def test_esp_cost_cube3d_invariance_origin():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
        get_random_esp_cost_cube3d_args()
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        shift = np.random.uniform(-3, 3, 3)
        tmp = coordinates+shift
        grid = UniformGrid(shift, grid_rvecs, shape, pbc)
        cost = ESPCost.from_grid_data(tmp, grid, vref, weights)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)


def test_esp_cost_cube0d_invariance_origin():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
        get_random_esp_cost_cube0d_args()
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        shift = np.random.uniform(-3, 3, 3)
        tmp = coordinates+shift
        grid = UniformGrid(shift, grid_rvecs, shape, pbc)
        cost = ESPCost.from_grid_data(tmp, grid, vref, weights)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)


def test_esp_cost_cube3d_invariance_rotation():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
        get_random_esp_cost_cube3d_args()
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        A = np.random.normal(0, 1, (3, 3))
        A = 0.5*(A+A.T)
        evals, evecs = np.linalg.eigh(A)
        del A
        del evals
        new_grid_rvecs = np.dot(grid_rvecs, evecs)
        new_coordinates = np.dot(coordinates-origin, evecs)+origin

        grid = UniformGrid(origin, new_grid_rvecs, shape, pbc)
        cost = ESPCost.from_grid_data(new_coordinates, grid, vref, weights)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)


def test_esp_cost_cube0d_invariance_rotation():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
        get_random_esp_cost_cube0d_args()
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        A = np.random.normal(0, 1, (3, 3))
        A = 0.5*(A+A.T)
        evals, evecs = np.linalg.eigh(A)
        del A
        del evals
        new_grid_rvecs = np.dot(grid_rvecs, evecs)
        new_coordinates = np.dot(coordinates-origin, evecs)+origin

        grid = UniformGrid(origin, new_grid_rvecs, shape, pbc)
        cost = ESPCost.from_grid_data(new_coordinates, grid, vref, weights)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)


def test_esp_cost_cube3d_invariance_images():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
        get_random_esp_cost_cube3d_args()
    grid = UniformGrid(origin, grid_rvecs, shape, pbc)
    # Generate costs with displaced origin
    costs = []
    rvecs = grid.get_cell().rvecs
    for i in xrange(10):
        new_coordinates = coordinates.copy()
        for j in xrange(len(coordinates)):
            new_coordinates[j] += np.dot(np.random.randint(-3, 4, 3), rvecs)

        cost = ESPCost.from_grid_data(new_coordinates, grid, vref, weights)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs)


def test_esp_cost_cube3d_invariance_rcut():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
        get_random_esp_cost_cube3d_args()
    grid = UniformGrid(origin, grid_rvecs, shape, pbc)
    # Generate costs with displaced origin
    costs = []
    for i in xrange(10):
        rcut = np.random.uniform(10, 30)
        alpha = 4.5/rcut
        gcut = 1.5*alpha
        cost = ESPCost.from_grid_data(coordinates, grid, vref, weights, rcut=rcut, alpha=alpha, gcut=gcut)
        costs.append(cost)
    # Compare the cost functions
    check_costs(costs, eps1=1e-8)


def test_esp_cost_cube3d_gradient():
    # Some parameters
    coordinates, numbers, origin, grid_rvecs, shape, pbc, vref, weights = \
        get_random_esp_cost_cube3d_args()
    grid = UniformGrid(origin, grid_rvecs, shape, pbc)
    cost = ESPCost.from_grid_data(coordinates, grid, vref, weights)

    x0 = np.random.uniform(-0.5, 0.5, len(numbers)+1)
    dxs = np.random.uniform(-1e-5, 1e-5, (100, len(numbers)+1))
    check_delta(cost.value, cost.gradient, x0, dxs)


def test_esp_cost_solve():
    A = np.random.uniform(-1, 1, (11, 11))
    A = np.dot(A, A.T)
    B = np.random.uniform(-1, 1, 11)
    C = np.random.uniform(-1, 1, ())
    cost = ESPCost(A, B, C, 10)

    # test without constraint
    x = cost.solve()
    assert abs(cost.gradient(x)).max() < 1e-8

    # test with constraint 0.0
    x = cost.solve(qtot=0.0)
    charges = x[:10]
    assert abs(charges.sum()) < 1e-8
    gradient = cost.gradient(x)
    assert abs(gradient[:10] - gradient[:10].mean()).max() < 1e-9
    assert abs(gradient[10:]).max() < 1e-9

    # test with constraint 1.0
    x = cost.solve(qtot=1.0)
    charges = x[:10]
    assert abs(charges.sum()-1) < 1e-9
    gradient = cost.gradient(x)
    assert abs(gradient[:10] - gradient[:10].mean()).max() < 1e-9
    assert abs(gradient[10:]).max() < 1e-9


def test_esp_cost_solve_regularized():
    A = np.random.uniform(-1, 1, (11, 11))
    A = np.dot(A, A.T)
    B = np.random.uniform(-1, 1, 11)
    C = np.random.uniform(-1, 1, ())
    cost = ESPCost(A, B, C, 10)

    # test without constraint
    ridge = 1e-6
    x = cost.solve(ridge=ridge)
    l = ridge*np.diag(cost._A)[:cost.natom].mean()
    ls = np.ones(cost._B.shape)*l
    ls[-1] = 0
    assert abs(cost.gradient(x) + 2*ls*x).max() < 1e-9

    # test with constraint 0.0
    x = cost.solve(qtot=0.0, ridge=ridge)
    charges = x[:10]
    assert abs(charges.sum()) < 1e-9
    gradient = cost.gradient(x) + 2*ls*x
    assert abs(gradient[:10] - gradient[:10].mean()).max() < 1e-10
    assert abs(gradient[10:]).max() < 1e-10

    # test with constraint 1.0
    x = cost.solve(qtot=1.0, ridge=ridge)
    charges = x[:10]
    assert abs(charges.sum()-1) < 1e-9
    gradient = cost.gradient(x) + 2*ls*x
    assert abs(gradient[:10] - gradient[:10].mean()).max() < 1e-10
    assert abs(gradient[10:]).max() < 1e-10


def test_compare_cubetools():
    # Load structure from cube file
    fn_cube = context.get_fn('test/jbw_coarse_aedens.cube')
    mol = Molecule.from_file(fn_cube)

    # Use a different grid
    origin = np.array([0.0, 0.0, 0.0])
    grid_rvecs = np.diag([0.183933, 0.187713, 0.190349])*3
    shape = np.array([18, 25, 27])
    pbc = np.array([1, 1, 1])
    ugrid = UniformGrid(origin, grid_rvecs, shape, pbc)

    # Generate weights
    weights = setup_weights(mol.coordinates, mol.numbers, ugrid,
        near={8: (1.8*angstrom, 0.5*angstrom),
              14: (1.8*angstrom, 0.5*angstrom)})
    weights /= weights.sum()

    # Random ref data
    esp = np.random.uniform(-1, 1, shape)

    # Cost function
    cost = ESPCost.from_grid_data(mol.coordinates, ugrid, esp, weights)

    # Sanity checks
    gvol = ugrid.get_grid_cell().volume
    assert abs(cost._A[-1, -1] - gvol) < 1e-8
    assert abs(cost._A.T - cost._A).max() < 1e-10

    # Compare numbers withreference values obtained with cfit2.cubetools
    assert abs(cost._A[0, 0] - 0.002211777956341868*gvol) < 1e-7
    assert abs(cost._A[7, 3] - 1.5124537688153498E-4*gvol) < 1e-8
    assert abs(cost._A[18, 2] - -0.028812329600683098*gvol) < 1e-6


def test_worst():
    cost = get_random_esp_cost_cube3d()
    N = cost.natom
    assert cost.worst() < cost._C
    assert cost.worst(1.0) < cost._C + cost._A[:N,:N].sum()/N**2 - 2*cost._B[:N].sum()/N
    assert cost.worst(-1.0) < cost._C + cost._A[:N,:N].sum()/N**2 + 2*cost._B[:N].sum()/N
    assert cost.worst() > 0.0
    assert cost.worst(1.0) > 0.0
    assert cost.worst(-1.0) > 0.0


def test_value_charges1():
    cost = get_random_esp_cost_cube3d()
    x = cost.solve()
    charges = x[:-1]
    assert abs(cost.value(x) - cost.value_charges(charges)) < 1e-10


def test_value_charges2():
    cost = get_random_esp_cost_cube3d()
    for i in xrange(100):
        x = np.random.normal(0, 1, len(cost._A))
        x -= x.sum()
        assert cost.value(x) >= cost.value_charges(x[:-1])


def test_consistent():
    # random system
    natom = 5
    coordinates = np.random.uniform(0, 10, (natom, 3))
    numbers = np.ones(natom, int)
    charges = np.random.normal(0, 1, natom)
    charges -= charges.mean()

    # random grid
    origin = np.random.uniform(-3, 3, 3)
    nrep = 10
    grid_rvecs = np.diag(np.random.uniform(8.0, 10.0, 3))/nrep
    grid_rvecs += np.random.uniform(-1.0, 1.0, (3, 3))/nrep
    shape = np.array([nrep, nrep, nrep])
    pbc = np.array([1, 1, 1])
    ugrid = UniformGrid(origin, grid_rvecs, shape, pbc)

    # compute the ESP
    rcut = 20.0
    alpha = 3.0/10.0
    gcut = 1.1*alpha
    esp = np.zeros(shape)
    compute_esp_grid_cube(ugrid, esp, coordinates, charges, rcut, alpha, gcut)

    # Set up weights
    weights = np.ones(shape)
    for i in xrange(natom):
        multiply_near_mask(coordinates[i], ugrid, 1.0, 0.5, weights)

    # Fit the charges and test
    cost = ESPCost.from_grid_data(coordinates, ugrid, esp, weights)
    x = cost.solve()
    assert cost.value_charges(charges) < 1e-7
    assert cost.value(x) < 1e-7
    assert abs(charges - x[:-1]).max() < 1e-4


def test_hdf5():
    with h5.File('horton.espfit.test.test_cost.test_hdf5.h5', driver='core', backing_store=False) as f:
        cost1 = get_random_esp_cost_cube3d()
        cost1.to_hdf5(f)
        cost2 = ESPCost.from_hdf5(f, None)
    assert abs(cost1._A - cost2._A).max() == 0
    assert abs(cost1._B - cost2._B).max() == 0
    assert cost1._C - cost2._C == 0
    assert cost1.natom == cost2.natom
