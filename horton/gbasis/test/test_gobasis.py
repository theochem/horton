# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


def test_shell_nbasis():
    assert get_shell_nbasis(-3) == 7
    assert get_shell_nbasis(-2) == 5
    assert get_shell_nbasis( 0) == 1
    assert get_shell_nbasis( 1) == 3
    assert get_shell_nbasis( 2) == 6
    assert get_shell_nbasis( 3) == 10
    try:
        get_shell_nbasis(-1)
        assert False
    except ValueError:
        pass


def test_gobasis_consistency():
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 0, 1, 1, 1, 1])
    nprims = np.array([2, 3, 3, 5, 5, 5, 7])
    shell_types = np.array([2, 1, 0, -2, 3, 0, 1])
    alphas = np.random.uniform(0, 1, nprims.sum())
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())

    gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    assert gobasis.nbasis == 29
    assert gobasis.max_shell_type == 3
    scales = gobasis.get_scales()
    assert abs(scales[0] - gob_normalization(alphas[0], np.array([2, 0, 0]))) < 1e-10

    shell_types = np.array([1, 1, 0, -2, -2, 0, 1])
    gobasis = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    assert gobasis.nbasis == 21
    assert gobasis.max_shell_type == 2

    # The center indexes in the shell_map are out of range.
    shell_map[0] = 2
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    shell_map[0] = 0

    # The size of the array shell_types does not match the sum of nprims.
    shell_types = np.array([1, 1])
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except TypeError:
        pass
    shell_types = np.array([1, 1, 0, -2, -2, 0, 1])

    # The elements of nprims should be at least 1.
    nprims[1] = 0
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    nprims[1] = 3

    # The size of the array alphas does not match the sum of nprims.
    alphas = np.random.uniform(-1, 1, 2)
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except TypeError:
        pass
    alphas = np.random.uniform(-1, 1, nprims.sum())

    # Encountered the nonexistent shell_type -1.
    shell_types[1] = -1
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    shell_types[1] = 1

    # The size of con_coeffs does not match nprims.
    con_coeffs = np.random.uniform(-1, 1, 3)
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except TypeError:
        pass
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())

    # Exceeding the maximym shell type (above):
    shell_types[0] = get_max_shell_type()+1
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    shell_types[0] = 2

    # Exceeding the maximym shell type (below):
    shell_types[0] = -get_max_shell_type()-1
    try:
        i2 = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except ValueError:
        pass
    shell_types[0] = 2


def test_load_basis():
    for go_basis_family in go_basis_families.itervalues():
        go_basis_family.load()


def test_grid_lih_321g_hf_some_points():
    ref = np.array([ # from cubegen
        [0.0, 0.0, 0.0, 0.037565082428],
        [0.1, 0.0, 0.0, 0.034775306876],
        [0.0, 0.1, 0.0, 0.034775306876],
        [0.0, 0.0, 1.0, 0.186234028507],
        [0.4, 0.2, 0.1, 0.018503681370],
    ])
    ref[:,:3] *= angstrom # convert from angstrom, cubegen docs are wrong. pfff.
    sys = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))

    # check for one point the compute_grid method
    output = np.zeros(sys.obasis.nbasis, float)
    point = np.array([0.0, 0.0, 1.0])*angstrom
    grid_fn = GB1GridDensityFn(sys.obasis.max_shell_type)
    sys.obasis.compute_grid1(output, point, grid_fn)
    # first basis function is contraction of three s-type gaussians
    assert sys.obasis.nprims[0] == 3
    scales = sys.obasis.get_scales()
    total = 0.0
    for i in xrange(3):
        alpha = sys.obasis.alphas[i]
        coeff = sys.obasis.con_coeffs[i]
        nrml = gob_normalization(alpha, np.zeros(3, int))
        # check scale
        assert abs(scales[i] - nrml) < 1e-10
        # check that we are on the first atom
        assert sys.obasis.shell_map[i] == 0
        dsq = np.linalg.norm(point - sys.coordinates[0])**2
        gauss = nrml*np.exp(-alpha*dsq)
        total += coeff*gauss
    assert abs(total - output[0]) < 1e-10

    # check density matrix value
    dm = sys.wfn.dm_full
    assert abs(dm._array[0,0] - 1.96589709) < 1e-7

    points = ref[:,:3].copy()
    rhos = sys.compute_grid_density(points)
    assert abs(rhos - ref[:,3]).max() < 1e-5


def test_gird_one_body_ne():
    sys = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 2e1, 100)
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False)
    dist0 = np.sqrt(((grid.points - sys.coordinates[0])**2).sum(axis=1))
    dist1 = np.sqrt(((grid.points - sys.coordinates[1])**2).sum(axis=1))
    pot = sys.numbers[0]/dist0 + sys.numbers[1]/dist1
    na_ana = sys.get_nuclear_attraction()
    na_grid = sys.lf.create_one_body()
    sys.compute_grid_density_fock(grid.points, grid.weights, pot, na_grid)
    assert abs(na_grid._array).max() > 8.0
    assert abs(na_ana._array-na_grid._array).max() < 2e-3
    # check symmetry
    na_grid.check_symmetry()
