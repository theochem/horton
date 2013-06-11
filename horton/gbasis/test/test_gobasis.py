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
    assert abs(scales[0] - gob_cart_normalization(alphas[0], np.array([2, 0, 0]))) < 1e-10

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


def test_grid_lih_321g_hf_density_some_points():
    ref = np.array([ # from cubegen
        [0.0, 0.0, 0.0, 0.037565082428],
        [0.1, 0.0, 0.0, 0.034775306876],
        [0.0, 0.1, 0.0, 0.034775306876],
        [0.0, 0.0, 1.0, 0.186234028507],
        [0.4, 0.2, 0.1, 0.018503681370],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))

    # check for one point the compute_grid_point1 method
    output = np.zeros(sys.obasis.nbasis, float)
    point = np.array([0.0, 0.0, 1.0])*angstrom
    grid_fn = GB1DMGridDensityFn(sys.obasis.max_shell_type)
    sys.obasis.compute_grid_point1(output, point, grid_fn)
    # first basis function is contraction of three s-type gaussians
    assert sys.obasis.nprims[0] == 3
    scales = sys.obasis.get_scales()
    total = 0.0
    for i in xrange(3):
        alpha = sys.obasis.alphas[i]
        coeff = sys.obasis.con_coeffs[i]
        nrml = gob_cart_normalization(alpha, np.zeros(3, int))
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


def test_grid_co_ccpv5z_cart_hf_density_some_points():
    ref = np.array([ # from cubegen
        [ 0.0,  0.0,  0.0,   4.54392441417],
        [ 0.1,  0.0,  0.0,   2.87874696902],
        [ 0.0,  0.1,  0.0,   2.90909931711],
        [ 0.0,  0.0,  1.0,   0.00563354926],
        [ 0.4,  0.2,  0.1,   0.15257439924],
        [-0.4,  0.2,  0.1,   0.14408104500],
        [ 0.4, -0.2,  0.1,   0.14627065655],
        [ 0.4,  0.2, -0.1,   0.11912840380],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'))
    points = ref[:,:3].copy()
    rhos = sys.compute_grid_density(points)
    assert abs(rhos - ref[:,3]).max() < 3e-3


def test_grid_co_ccpv5z_pure_hf_density_some_points():
    ref = np.array([ # from cubegen
        [ 0.0,  0.0,  0.0,   4.54338939220],
        [ 0.1,  0.0,  0.0,   2.87742753163],
        [ 0.0,  0.1,  0.0,   2.90860415538],
        [ 0.0,  0.0,  1.0,   0.00285462032],
        [ 0.4,  0.2,  0.1,   0.15399703660],
        [-0.4,  0.2,  0.1,   0.14425254494],
        [ 0.4, -0.2,  0.1,   0.14409038614],
        [ 0.4,  0.2, -0.1,   0.11750780363],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'))
    points = ref[:,:3].copy()
    rhos = sys.compute_grid_density(points)
    assert abs(rhos - ref[:,3]).max() < 3e-3


def test_grid_lih_321g_hf_gradient_some_points():
    ref = np.array([ # from cubegen
        [0.0, 0.0, 0.0,  0.000000000000,  0.000000000000,  0.179349665782],
        [0.1, 0.0, 0.0, -0.028292898754,  0.000000000000,  0.164582727812],
        [0.0, 0.1, 0.0,  0.000000000000, -0.028292898754,  0.164582727812],
        [0.0, 0.0, 1.0,  0.000000000000,  0.000000000000, -0.929962409854],
        [0.4, 0.2, 0.1, -0.057943497876, -0.028971748938,  0.069569174116],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
    points = ref[:,:3].copy()
    gradients = sys.compute_grid_gradient(points)
    assert abs(gradients - ref[:,3:]).max() < 1e-6


def test_grid_co_ccpv5z_cart_hf_gradient_some_points():
    ref = np.array([ # from cubegen
        [ 0.0,  0.0,  0.0,  -0.26805895992,  -0.03725931097,  26.06939895580],
        [ 0.1,  0.0,  0.0, -11.66097634913,  -0.02427222636,  11.49946087301],
        [ 0.0,  0.1,  0.0,  -0.18730587145, -11.60371334591,  11.60046471817],
        [ 0.0,  0.0,  1.0,   0.00350647376,  -0.00151630329,  -0.00944412097],
        [ 0.4,  0.2,  0.1,  -0.46814335442,  -0.28380627268,  -0.02592227656],
        [-0.4,  0.2,  0.1,   0.63742782898,  -0.32989678808,   0.00444361306],
        [ 0.4, -0.2,  0.1,  -0.50464249640,   0.29978538874,  -0.01244489023],
        [ 0.4,  0.2, -0.1,  -0.21837773815,  -0.16855926400,   0.15518115326],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'))
    points = ref[:,:3].copy()
    gradrhos = sys.compute_grid_gradient(points)
    assert abs(gradrhos - ref[:,3:]).max() < 1e-2 # cubegen output somehow not reliable


def test_grid_co_ccpv5z_pure_hf_gradient_some_points():
    ref = np.array([ # from cubegen
        [ 0.0,  0.0,  0.0,  -0.27796827654,  -0.03971005800,  26.06788123216],
        [ 0.1,  0.0,  0.0, -11.65999871789,  -0.02706024561,  11.49763108605],
        [ 0.0,  0.1,  0.0,  -0.19499030621, -11.60235682832,  11.60235521243],
        [ 0.0,  0.0,  1.0,   0.00184843964,   0.00026806115,  -0.01003272687],
        [ 0.4,  0.2,  0.1,  -0.46500454519,  -0.27516942731,  -0.01707049479],
        [-0.4,  0.2,  0.1,   0.63911725484,  -0.32989616481,   0.00229353087],
        [ 0.4, -0.2,  0.1,  -0.51099806603,   0.29961935521,  -0.00979594206],
        [ 0.4,  0.2, -0.1,  -0.21849813344,  -0.16098019809,   0.16093849962],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'))
    points = ref[:,:3].copy()
    gradrhos = sys.compute_grid_gradient(points)
    assert abs(gradrhos - ref[:,3:]).max() < 1e-4


def test_grid_lih_321g_hf_esp_some_points():
    ref = np.array([ # from cubegen
        [0.0, 0.0, 0.0, 0.906151727538],
        [0.1, 0.0, 0.0, 0.891755005233],
        [0.0, 0.1, 0.0, 0.891755005233],
        [0.0, 0.0, 1.0, 1.422294470114],
        [0.4, 0.2, 0.1, 0.796490099689],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))

    # check for one point the compute_grid_point2 method
    if False: # TODO Write python interface to GB2GridHartreeFn
    #for row in ref:
        point = row[:3]
        grid_fn = GB2GridHartreeFn(sys.obasis.max_shell_type)
        esp = sys.obasis.compute_grid_point2(sys.wfn.dm_full, point, grid_fn)
        hartree = -esp
        for n, pos in zip(sys.numbers, sys.coordinates):
            hartree -= n/np.linalg.norm(pos - point)
        assert abs(row[3] - hartree) < 1e-8

    points = ref[:,:3].copy()
    esps = sys.compute_grid_esp(points)
    assert abs(esps - ref[:,3]).max() < 1e-8


def test_grid_co_ccpv5z_cart_hf_esp_some_points():
    ref = np.array([ # from cubegen
        [ 0.0,  0.0,  0.0,  10.69443507172],
        [ 0.1,  0.0,  0.0,   6.43122889229],
        [ 0.0,  0.1,  0.0,   6.43406765938],
        [ 0.0,  0.0,  1.0,   0.27023448629],
        [ 0.4,  0.2,  0.1,   0.82646540602],
        [-0.4,  0.2,  0.1,   0.93595072191],
        [ 0.4, -0.2,  0.1,   0.83432301119],
        [ 0.4,  0.2, -0.1,   0.68524674809],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'))
    points = ref[:,:3].copy()
    esps = sys.compute_grid_esp(points)
    assert abs(esps - ref[:,3]).max() < 1e-3 # cubegen output somehow not reliable


def test_grid_co_ccpv5z_pure_hf_esp_some_points():
    ref = np.array([ # from cubegen
        [ 0.0,  0.0,  0.0,  10.69443507172],
        [ 0.1,  0.0,  0.0,   6.43122889229],
        [ 0.0,  0.1,  0.0,   6.43406765938],
        [ 0.0,  0.0,  1.0,   0.27023448629],
        [ 0.4,  0.2,  0.1,   0.82646540602],
        [-0.4,  0.2,  0.1,   0.93595072191],
        [ 0.4, -0.2,  0.1,   0.83432301119],
        [ 0.4,  0.2, -0.1,   0.68524674809],
    ])
    ref[:,:3] *= angstrom
    sys = System.from_file(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'))
    points = ref[:,:3].copy()
    esps = sys.compute_grid_esp(points)
    assert abs(esps - ref[:,3]).max() < 1e-5


def test_grid_one_body_ne():
    sys = System.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
    rtf = ExpRTransform(1e-3, 2e1, 100)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False)
    dist0 = np.sqrt(((grid.points - sys.coordinates[0])**2).sum(axis=1))
    dist1 = np.sqrt(((grid.points - sys.coordinates[1])**2).sum(axis=1))
    pot = -sys.numbers[0]/dist0 - sys.numbers[1]/dist1
    na_ana = sys.get_nuclear_attraction()
    na_grid = sys.lf.create_one_body()
    sys.compute_grid_density_fock(grid.points, grid.weights, pot, na_grid)
    assert abs(na_grid._array).max() > 8.0
    assert abs(na_ana._array-na_grid._array).max() < 2e-3
    # check symmetry
    na_grid.check_symmetry()


def test_gob_normalization():
    assert abs(gob_pure_normalization(0.09515, 0) - 0.122100288) < 1e-5
    assert abs(gob_pure_normalization(0.1687144, 1) - 0.154127551) < 1e-5
    assert abs(gob_cart_normalization(0.344, np.array([1,1,0])) - 0.440501466) < 1e-8
    assert abs(gob_cart_normalization(0.246, np.array([1,1,1])) - 0.242998767) < 1e-8
    assert abs(gob_cart_normalization(0.238, np.array([2,1,1])) - 0.127073818) < 1e-8
    assert abs(gob_pure_normalization(0.3, 0) - gob_cart_normalization(0.3, np.array([0, 0, 0]))) < 1e-10
    assert abs(gob_pure_normalization(0.7, 0) - gob_cart_normalization(0.7, np.array([0, 0, 0]))) < 1e-10
    assert abs(gob_pure_normalization(1.9, 0) - gob_cart_normalization(1.9, np.array([0, 0, 0]))) < 1e-10
    assert abs(gob_pure_normalization(0.3, 1) - gob_cart_normalization(0.3, np.array([1, 0, 0]))) < 1e-10
    assert abs(gob_pure_normalization(0.7, 1) - gob_cart_normalization(0.7, np.array([0, 1, 0]))) < 1e-10
    assert abs(gob_pure_normalization(1.9, 1) - gob_cart_normalization(1.9, np.array([0, 0, 1]))) < 1e-10


def test_cart_pure_switch():
    sys = System.from_file(context.get_fn('test/water.xyz'), obasis='aug-cc-pvdz')
    assert sys.obasis.nbasis == 41
    sys = System.from_file(context.get_fn('test/water.xyz'), obasis=GOBasisDesc('aug-cc-pvdz', pure=False))
    assert sys.obasis.nbasis == 43


def get_olp(ob):
    lf = DenseLinalgFactory(ob.nbasis)
    olp = lf.create_one_body()
    ob.compute_overlap(olp)
    return olp._array

def test_concatenate1():
    sys = System.from_file(context.get_fn('test/water.xyz'), obasis='3-21g')
    ob = GOBasis.concatenate(sys.obasis, sys.obasis)
    assert ob.ncenter == 3*2
    assert ob.nbasis == 13*2
    a = get_olp(ob)
    assert abs(a[:13,:13] - a[:13,13:]).max() < 1e-15
    assert (a[:13,:13] == a[13:,13:]).all()
    assert abs(a[:13,:13] - a[13:,:13]).max() < 1e-15


def test_concatenate2():
    sys1 = System.from_file(context.get_fn('test/water.xyz'), obasis='3-21g')
    sys2 = System.from_file(context.get_fn('test/water.xyz'), obasis='sto-3g')
    obasis = GOBasis.concatenate(sys1.obasis, sys2.obasis)
    assert obasis.ncenter == 3*2
    assert obasis.nbasis == sys1.obasis.nbasis+sys2.obasis.nbasis

    a = get_olp(obasis)
    a11 = get_olp(sys1.obasis)
    a22 = get_olp(sys2.obasis)
    N = sys1.obasis.nbasis
    assert (a[:N,:N] == a11).all()
    assert (a[N:,N:] == a22).all()


def test_abstract():
    try:
        centers = np.zeros((1,3), float)
        shell_map = np.zeros(2, int)
        nprims = np.array([1, 2])
        shell_types = np.array([0, 1])
        alphas = np.array([1.0, 1.1, 1.2])
        con_coeffs = np.array([0.1, 0.2, 0.3])
        from horton.gbasis.cext import GBasis
        gb = GBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
        assert False
    except NotImplementedError:
        pass
