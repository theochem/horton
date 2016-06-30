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

from horton.grid.test.common import get_cosine_spline
from horton.test.common import get_random_cell, numpy_seed


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


def test_grid_integrate_segments():
    npoint = 10
    segments = np.array([2, 5, 3])
    grid = IntGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    pot = np.random.normal(0, 1, npoint)
    dens = np.random.normal(0, 1, npoint)

    # three
    ints = grid.integrate(pot, dens, segments=segments)
    assert ints.shape == (3,)
    product = grid.weights*pot*dens
    assert abs(ints[0] - product[:2].sum()) < 1e-10
    assert abs(ints[1] - product[2:7].sum()) < 1e-10
    assert abs(ints[2] - product[7:].sum()) < 1e-10

    # two
    ints = grid.integrate(pot, segments=segments)
    assert ints.shape == (3,)
    product = grid.weights*pot
    assert abs(ints[0] - product[:2].sum()) < 1e-10
    assert abs(ints[1] - product[2:7].sum()) < 1e-10
    assert abs(ints[2] - product[7:].sum()) < 1e-10

    # one
    ints = grid.integrate(segments=segments)
    assert ints.shape == (3,)
    assert abs(ints[0] - grid.weights[:2].sum()) < 1e-10
    assert abs(ints[1] - grid.weights[2:7].sum()) < 1e-10
    assert abs(ints[2] - grid.weights[7:].sum()) < 1e-10


def test_grid_integrate_cartesian_moments():
    npoint = 10
    grid = IntGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    dens = np.random.normal(0, 1, npoint)

    center = np.random.normal(0, 1, 3)
    x = grid.points[:,0]-center[0]
    y = grid.points[:,1]-center[1]
    z = grid.points[:,2]-center[2]

    ints = grid.integrate(dens, center=center, lmax=2, mtype=1)
    assert ints.shape == (10,)
    assert abs(ints[0] - (grid.weights*dens).sum()) < 1e-10
    assert abs(ints[1] - (grid.weights*dens*x).sum()) < 1e-10
    assert abs(ints[2] - (grid.weights*dens*y).sum()) < 1e-10
    assert abs(ints[3] - (grid.weights*dens*z).sum()) < 1e-10
    assert abs(ints[4] - (grid.weights*dens*x*x).sum()) < 1e-10
    assert abs(ints[5] - (grid.weights*dens*x*y).sum()) < 1e-10
    assert abs(ints[9] - (grid.weights*dens*z*z).sum()) < 1e-10


def test_grid_integrate_cartesian_moments_segments():
    npoint = 10
    segments = np.array([2, 5, 3])
    grid = IntGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    dens = np.random.normal(0, 1, npoint)

    center = np.random.normal(0, 1, 3)
    x = grid.points[:,0]-center[0]
    y = grid.points[:,1]-center[1]
    z = grid.points[:,2]-center[2]

    ints = grid.integrate(dens, center=center, lmax=2, mtype=1, segments=segments)
    assert ints.shape == (3, 10)
    for i, begin, end in (0, 0, 2), (1, 2, 7), (2, 7, 10):
        assert abs(ints[i, 0] - (grid.weights*dens)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 1] - (grid.weights*dens*x)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 2] - (grid.weights*dens*y)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 3] - (grid.weights*dens*z)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 4] - (grid.weights*dens*x*x)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 5] - (grid.weights*dens*x*y)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 9] - (grid.weights*dens*z*z)[begin:end].sum()) < 1e-10


def test_grid_integrate_pure_moments():
    npoint = 10
    grid = IntGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    dens = np.random.normal(0, 1, npoint)

    center = np.random.normal(0, 1, 3)
    x = grid.points[:,0]-center[0]
    y = grid.points[:,1]-center[1]
    z = grid.points[:,2]-center[2]
    r2 = x*x+y*y+z*z

    ints = grid.integrate(dens, center=center, lmax=2, mtype=2)
    assert ints.shape == (9,)
    assert abs(ints[0] - (grid.weights*dens).sum()) < 1e-10
    assert abs(ints[1] - (grid.weights*dens*z).sum()) < 1e-10
    assert abs(ints[2] - (grid.weights*dens*x).sum()) < 1e-10
    assert abs(ints[3] - (grid.weights*dens*y).sum()) < 1e-10
    assert abs(ints[4] - (grid.weights*dens*(1.5*z**2-0.5*r2)).sum()) < 1e-10
    assert abs(ints[5] - (grid.weights*dens*(3.0**0.5*x*z)).sum()) < 1e-10
    assert abs(ints[8] - (grid.weights*dens*(3.0**0.5*x*y)).sum()) < 1e-10


def test_grid_integrate_pure_moments_segments():
    npoint = 10
    segments = np.array([2, 5, 3])
    grid = IntGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    dens = np.random.normal(0, 1, npoint)

    center = np.random.normal(0, 1, 3)
    x = grid.points[:,0]-center[0]
    y = grid.points[:,1]-center[1]
    z = grid.points[:,2]-center[2]
    r2 = x*x+y*y+z*z

    ints = grid.integrate(dens, center=center, lmax=2, mtype=2, segments=segments)
    assert ints.shape == (3, 9)
    for i, begin, end in (0, 0, 2), (1, 2, 7), (2, 7, 10):
        assert abs(ints[i, 0] - (grid.weights*dens)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 1] - (grid.weights*dens*z)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 2] - (grid.weights*dens*x)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 3] - (grid.weights*dens*y)[begin:end].sum()) < 1e-10
        assert abs(ints[i, 4] - (grid.weights*dens*(1.5*z**2-0.5*r2))[begin:end].sum()) < 1e-10
        assert abs(ints[i, 5] - (grid.weights*dens*(3.0**0.5*x*z))[begin:end].sum()) < 1e-10
        assert abs(ints[i, 8] - (grid.weights*dens*(3.0**0.5*x*y))[begin:end].sum()) < 1e-10


def test_grid_integrate_radial_moments():
    npoint = 10
    grid = IntGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    dens = np.random.normal(0, 1, npoint)

    center = np.random.normal(0, 1, 3)
    x = grid.points[:,0]-center[0]
    y = grid.points[:,1]-center[1]
    z = grid.points[:,2]-center[2]
    r = np.sqrt(x*x+y*y+z*z)

    ints = grid.integrate(dens, center=center, lmax=2, mtype=3)
    ints.shape == (3,)
    assert abs(ints[0] - (grid.weights*dens).sum()) < 1e-10
    assert abs(ints[1] - (grid.weights*dens*r).sum()) < 1e-10
    assert abs(ints[2] - (grid.weights*dens*r*r).sum()) < 1e-10


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

        distances = np.sqrt(((center - g.points)**2).sum(axis=1))
        output2 = cs(distances)

        assert abs(output1 - output2).max() < 1e-10


def test_eval_spline_grid_add_random():
    npoint = 10
    cs = get_cosine_spline()

    for irep in xrange(10):
        with numpy_seed(irep):
            cell = get_random_cell(1.0, irep % 4)
            points = np.random.normal(-2, 3, (npoint,3))
            g = IntGrid(points, np.random.normal(0, 1.0, npoint))
            center1 = np.random.uniform(-2, 2, 3)
            center2 = np.random.uniform(-2, 2, 3)

        output1 = np.zeros(npoint)
        g.eval_spline(cs, center1, output1, cell)

        output2 = np.zeros(npoint)
        g.eval_spline(cs, center2, output2, cell)

        output3 = np.zeros(npoint)
        g.eval_spline(cs, center1, output3, cell)
        g.eval_spline(cs, center2, output3, cell)

        assert abs(output1 + output2 - output3).max() < 1e-10


def test_density_decomposition_n2():
    # compute reference density and becke_weights for the first atom
    mol = IOData.from_file(context.get_fn('test/n2_hfs_sto3g.fchk'))
    molgrid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'veryfine', random_rotate=False, mode='only')
    atgrid = molgrid.subgrids[0]
    dm_full = mol.get_dm_full()
    rho = mol.obasis.compute_grid_density_dm(dm_full, atgrid.points)
    becke_weights = np.ones(atgrid.size)
    becke_helper_atom(atgrid.points, becke_weights, np.ones(2), mol.coordinates, 0, 3)
    rho_ref = rho*becke_weights

    # generate real regular solid harmonics
    lmaxmax = 10
    delta = atgrid.points - atgrid.center
    x, y, z = delta.T
    r = np.sqrt(x*x + y*y + z*z)
    checks = {}

    work = np.zeros((atgrid.size, (lmaxmax+1)**2), float)
    work[:,0] = z/r
    work[:,1] = x/r
    work[:,2] = y/r
    fill_pure_polynomials(work, lmaxmax)

    # create decompositions of the Becke AIM density of the first atom
    splines = atgrid.get_spherical_decomposition(rho_ref[:atgrid.size], lmax=lmaxmax)

    # manual reconstruction for differen lmax
    tmp = np.zeros(atgrid.size)
    atgrid.eval_spline(splines[0], mol.coordinates[0], tmp)
    tmp /= np.sqrt(4*np.pi)
    checks[0] = tmp.copy()
    counter = 0
    for l in xrange(1, lmaxmax+1):
        for m in xrange(-l, l+1):
            output = np.zeros(atgrid.size)
            atgrid.eval_spline(splines[counter+1], mol.coordinates[0], output)
            output *= work[:,counter]/np.sqrt(4*np.pi)*np.sqrt(2*l+1)
            tmp += output
            counter += 1
        checks[l] = tmp.copy()
    assert counter == (lmaxmax+1)**2-1

    # reconstruct the atom at different levels of accuarcy
    last_error = None
    for lmax in xrange(0, lmaxmax+1):
        output = np.zeros(atgrid.size)
        atgrid.eval_decomposition(splines[:(lmax+1)**2], mol.coordinates[0], output)
        if lmax in checks:
            error = np.sqrt(atgrid.integrate((output - checks[lmax])**2))
            assert abs(error) < 1e-12
        error = np.sqrt(atgrid.integrate((output - rho_ref)**2))
        if last_error is not None:
            assert last_error > error
        last_error = error
    assert error < 1e-3

    # test sensibility of results at the positions of the nuclei. Here, the potentials
    # should be finite and non-zero.
    nucgrid = IntGrid(mol.coordinates, np.ones(mol.natom))
    # Contribution from s should be non-zero:
    tmp_s = np.zeros(nucgrid.size)
    nucgrid.eval_spline(splines[0], mol.coordinates[0], tmp_s)
    assert tmp_s[0] != 0.0
    np.testing.assert_almost_equal(tmp_s[0], splines[0](np.array([0.0]))[0])
    assert np.isfinite(tmp_s).all()
    # Contribution from p and higher should be zero.
    # This tested by computing all and comparing to contribution from s
    # (with proper scale factor).
    for lmax in xrange(0, lmaxmax+1):
        tmp = np.zeros(nucgrid.size)
        nucgrid.eval_decomposition(splines[:(lmax+1)**2], mol.coordinates[0], tmp)
        np.testing.assert_almost_equal(tmp[0], tmp_s[0]/np.sqrt(4*np.pi))
        assert np.isfinite(tmp).all()
