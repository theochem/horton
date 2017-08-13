# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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
import os
from nose.tools import assert_raises
from nose.plugins.attrib import attr

from .lightgrid import generate_molecular_grid
from .common import load_obasis, load_mdata, load_dm, load_orbsa_coeffs
from .. import (GOBasisDesc, GOBasis, GB1DMGridDensityFn, get_shell_nbasis, get_gobasis,
                gob_pure_normalization, gob_cart_normalization, go_basis_families,
                get_max_shell_type)


angstrom = 1.0e-10 / 0.5291772083e-10


def test_shell_nbasis():
    assert get_shell_nbasis(-3) == 7
    assert get_shell_nbasis(-2) == 5
    assert get_shell_nbasis(0) == 1
    assert get_shell_nbasis(1) == 3
    assert get_shell_nbasis(2) == 6
    assert get_shell_nbasis(3) == 10
    with assert_raises(ValueError):
        get_shell_nbasis(-1)


def test_gobasis_consistency():
    centers = np.random.uniform(-1, 1, (2, 3))
    shell_map = np.array([0, 0, 0, 1, 1, 1, 1])
    nprims = np.array([2, 3, 3, 5, 5, 5, 7])
    shell_types = np.array([2, 1, 0, -2, 3, 0, 1])
    alphas = np.random.uniform(0, 1, nprims.sum())
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())

    gb = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    assert gb.nbasis == 29
    assert gb.max_shell_type == 3
    scales = gb.get_scales()
    assert abs(scales[0] - gob_cart_normalization(alphas[0], np.array([2, 0, 0]))) < 1e-10
    assert (gb.basis_offsets == np.array([0, 6, 9, 10, 15, 25, 26])).all()
    assert (gb.shell_lookup == np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 3,
                                         3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4,
                                         4, 4, 4, 5, 6, 6, 6])).all()

    shell_types = np.array([1, 1, 0, -2, -2, 0, 1])
    gb = GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    assert gb.nbasis == 21
    assert gb.max_shell_type == 2

    # The center indexes in the shell_map are out of range.
    shell_map[0] = 2
    with assert_raises(ValueError):
        GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    shell_map[0] = 0

    # The size of the array shell_types does not match the sum of nprims.
    shell_types = np.array([1, 1])
    with assert_raises(TypeError):
        GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    shell_types = np.array([1, 1, 0, -2, -2, 0, 1])

    # The elements of nprims should be at least 1.
    nprims[1] = 0
    with assert_raises(ValueError):
        GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    nprims[1] = 3

    # The size of the array alphas does not match the sum of nprims.
    alphas = np.random.uniform(-1, 1, 2)
    with assert_raises(TypeError):
        GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    alphas = np.random.uniform(-1, 1, nprims.sum())

    # Encountered the nonexistent shell_type -1.
    shell_types[1] = -1
    with assert_raises(ValueError):
        GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    shell_types[1] = 1

    # The size of con_coeffs does not match nprims.
    con_coeffs = np.random.uniform(-1, 1, 3)
    with assert_raises(TypeError):
        GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    con_coeffs = np.random.uniform(-1, 1, nprims.sum())

    # Exceeding the maximym shell type (above):
    shell_types[0] = get_max_shell_type() + 1
    with assert_raises(ValueError):
        GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    shell_types[0] = 2

    # Exceeding the maximym shell type (below):
    shell_types[0] = -get_max_shell_type() - 1
    with assert_raises(ValueError):
        GOBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)
    shell_types[0] = 2


def test_load_basis():
    for go_basis_family in go_basis_families.values():
        assert os.path.basename(go_basis_family.filename).islower()
        go_basis_family.load()


def test_grid_lih_321g_hf_density_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, 0.037565082428],
        [0.1, 0.0, 0.0, 0.034775306876],
        [0.0, 0.1, 0.0, 0.034775306876],
        [0.0, 0.0, 1.0, 0.186234028507],
        [0.4, 0.2, 0.1, 0.018503681370],
    ])
    ref[:, :3] *= angstrom
    fn = 'li_h_3_21G_hf_g09_fchk'
    obasis = load_obasis(fn)
    mol = load_mdata(fn)

    # check for one point the compute_grid_point1 method
    output = np.zeros(obasis.nbasis, float)
    point = np.array([0.0, 0.0, 1.0]) * angstrom
    grid_fn = GB1DMGridDensityFn(obasis.max_shell_type)
    obasis.compute_grid_point1(output, point, grid_fn)
    # first basis function is contraction of three s-type gaussians
    assert obasis.nprims[0] == 3
    scales = obasis.get_scales()
    total = 0.0
    for i in range(3):
        alpha = obasis.alphas[i]
        coeff = obasis.con_coeffs[i]
        nrml = gob_cart_normalization(alpha, np.zeros(3, int))
        # check scale
        assert abs(scales[i] - nrml) < 1e-10
        # check that we are on the first atom
        assert obasis.shell_map[i] == 0
        dsq = np.linalg.norm(point - mol['coordinates'][0]) ** 2
        gauss = nrml * np.exp(-alpha * dsq)
        total += coeff * gauss
    assert abs(total - output[0]) < 1e-10

    # check density matrix value
    dm_full = load_dm(fn)
    assert abs(dm_full[0, 0] - 1.96589709) < 1e-7

    points = ref[:, :3].copy()
    rhos = obasis.compute_grid_density_dm(dm_full, points)
    assert abs(rhos - ref[:, 3]).max() < 1e-5


def check_grid_rho(fn, ref, eps):
    points = ref[:, :3].copy()
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)
    rhos = obasis.compute_grid_density_dm(dm_full, points)
    assert abs(rhos - ref[:, 3]).max() < eps


def test_grid_co_ccpv5z_cart_hf_density_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, 4.54392441417],
        [0.1, 0.0, 0.0, 2.87874696902],
        [0.0, 0.1, 0.0, 2.90909931711],
        [0.0, 0.0, 1.0, 0.00563354926],
        [0.4, 0.2, 0.1, 0.15257439924],
        [-0.4, 0.2, 0.1, 0.14408104500],
        [0.4, -0.2, 0.1, 0.14627065655],
        [0.4, 0.2, -0.1, 0.11912840380],
    ])
    ref[:, :3] *= angstrom
    check_grid_rho('co_ccpv5z_cart_hf_g03_fchk', ref, 3e-3)


def test_grid_co_ccpv5z_pure_hf_density_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, 4.54338939220],
        [0.1, 0.0, 0.0, 2.87742753163],
        [0.0, 0.1, 0.0, 2.90860415538],
        [0.0, 0.0, 1.0, 0.00285462032],
        [0.4, 0.2, 0.1, 0.15399703660],
        [-0.4, 0.2, 0.1, 0.14425254494],
        [0.4, -0.2, 0.1, 0.14409038614],
        [0.4, 0.2, -0.1, 0.11750780363],
    ])
    ref[:, :3] *= angstrom
    check_grid_rho('co_ccpv5z_pure_hf_g03_fchk', ref, 3e-3)


def check_grid_gradient(fn, ref, eps):
    points = ref[:, :3].copy()
    obasis = load_obasis(fn)
    dm_full = load_dm(fn)
    gradients = obasis.compute_grid_gradient_dm(dm_full, points)
    assert abs(gradients - ref[:, 3:]).max() < eps


def test_grid_lih_321g_hf_gradient_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, 0.000000000000, 0.000000000000, 0.179349665782],
        [0.1, 0.0, 0.0, -0.028292898754, 0.000000000000, 0.164582727812],
        [0.0, 0.1, 0.0, 0.000000000000, -0.028292898754, 0.164582727812],
        [0.0, 0.0, 1.0, 0.000000000000, 0.000000000000, -0.929962409854],
        [0.4, 0.2, 0.1, -0.057943497876, -0.028971748938, 0.069569174116],
    ])
    ref[:, :3] *= angstrom
    check_grid_gradient('li_h_3_21G_hf_g09_fchk', ref, 1e-6)


def test_grid_lih_321g_hf_orbital_gradient_some_points():
    points = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [0.0, 0.1, 0.0],
        [0.0, 0.0, 1.0],
        [0.4, 0.2, 0.1]
    ])
    ref = np.array([  # calculated using finite difference
        [[0.000000000000, 0.000000000000, 0.335358022388],
         [0.000000000000, 0.000000000000, -0.039601719339],
         [0.000000000000, 0.000000000000, -0.065908783181],
         [0.066778303100, 0.000000000000, 0.000000000000],
         [0.000000000000, 0.066778303100, 0.000000000000],
         [0.000000000000, 0.000000000000, 0.035174712603],
         [0.000000000000, 0.000000000000, -0.086099568163],
         [0.121755587880, 0.000000000000, 0.000000000000],
         [0.000000000000, 0.121755587880, 0.000000000000],
         [0.000000000000, 0.000000000000, -0.068136399228],
         [0.000000000000, 0.000000000000, -0.027336484452]],
        [[-0.029716570085, 0.000000000000, 0.331085502737],
         [0.000124168300, 0.000000000000, -0.039374354768],
         [0.008077765337, 0.000000000000, -0.064966718776],
         [0.066265798547, 0.000000000000, 0.003804006584],
         [0.000000000000, 0.066607211612, 0.000000000000],
         [0.000377488384, 0.000000000000, 0.034880919059],
         [-0.000705284103, 0.000000000000, -0.085471067710],
         [0.120527897959, 0.000000000000, 0.009112195828],
         [0.000000000000, 0.121345725857, 0.000000000000],
         [0.013285776833, 0.000000000000, -0.067340615153],
         [-0.000503096838, 0.000000000000, -0.027343807139]],
        [[0.000000000000, -0.029716570085, 0.331085502737],
         [0.000000000000, 0.000124168300, -0.039374354768],
         [0.000000000000, 0.008077765337, -0.064966718776],
         [0.066607211612, 0.000000000000, 0.000000000000],
         [0.000000000000, 0.066265798547, 0.003804006584],
         [0.000000000000, 0.000377488384, 0.034880919059],
         [0.000000000000, -0.000705284103, -0.085471067710],
         [0.121345725857, 0.000000000000, 0.000000000000],
         [0.000000000000, 0.120527897959, 0.009112195828],
         [0.000000000000, 0.013285776833, -0.067340615153],
         [0.000000000000, -0.000503096838, -0.027343807139]],
        [[0.000000000000, 0.000000000000, 5.083843181919],
         [0.000000000000, 0.000000000000, -0.220159704611],
         [0.000000000000, 0.000000000000, -0.941775028074],
         [0.095151155766, 0.000000000000, 0.000000000000],
         [0.000000000000, 0.095151155766, 0.000000000000],
         [0.000000000000, 0.000000000000, 0.427204143785],
         [0.000000000000, 0.000000000000, -0.520952693280],
         [0.190225965755, 0.000000000000, 0.000000000000],
         [0.000000000000, 0.190225965755, 0.000000000000],
         [0.000000000000, 0.000000000000, -0.259089767114],
         [0.000000000000, 0.000000000000, 0.108036467562]],
        [[-0.122409907465, -0.061204953733, 0.310330235810],
         [0.001247852781, 0.000623926390, -0.037125981310],
         [0.032675748998, 0.016337874499, -0.058852274411],
         [0.061479984281, -0.002759557877, 0.013993646107],
         [-0.002759557877, 0.065619321096, 0.006996823053],
         [0.000434364406, 0.000217182203, 0.042093481225],
         [-0.000016445047, -0.000008222524, -0.094292047809],
         [0.109060186207, -0.006612206455, 0.033530326679],
         [-0.006612206454, 0.118978495887, 0.016765163340],
         [0.051799568647, 0.025899784323, -0.049145171706],
         [-0.001806896115, -0.000903448057, -0.019856823712]]
    ])
    fn = 'li_h_3_21G_hf_g09_fchk'
    obasis = load_obasis(fn)
    orb_alpha = load_orbsa_coeffs(fn)
    orbs = np.arange(obasis.nbasis)
    test = np.array(obasis.compute_grid_orb_gradient_exp(orb_alpha, points, orbs))

    np.testing.assert_almost_equal(test, ref, decimal=7)


def test_grid_co_ccpv5z_cart_hf_gradient_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, -0.26805895992, -0.03725931097, 26.06939895580],
        [0.1, 0.0, 0.0, -11.66097634913, -0.02427222636, 11.49946087301],
        [0.0, 0.1, 0.0, -0.18730587145, -11.60371334591, 11.60046471817],
        [0.0, 0.0, 1.0, 0.00350647376, -0.00151630329, -0.00944412097],
        [0.4, 0.2, 0.1, -0.46814335442, -0.28380627268, -0.02592227656],
        [-0.4, 0.2, 0.1, 0.63742782898, -0.32989678808, 0.00444361306],
        [0.4, -0.2, 0.1, -0.50464249640, 0.29978538874, -0.01244489023],
        [0.4, 0.2, -0.1, -0.21837773815, -0.16855926400, 0.15518115326],
    ])
    ref[:, :3] *= angstrom
    # cubegen output somehow not reliable?
    check_grid_gradient('co_ccpv5z_cart_hf_g03_fchk', ref, 1e-2)


def test_grid_co_ccpv5z_pure_hf_gradient_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, -0.27796827654, -0.03971005800, 26.06788123216],
        [0.1, 0.0, 0.0, -11.65999871789, -0.02706024561, 11.49763108605],
        [0.0, 0.1, 0.0, -0.19499030621, -11.60235682832, 11.60235521243],
        [0.0, 0.0, 1.0, 0.00184843964, 0.00026806115, -0.01003272687],
        [0.4, 0.2, 0.1, -0.46500454519, -0.27516942731, -0.01707049479],
        [-0.4, 0.2, 0.1, 0.63911725484, -0.32989616481, 0.00229353087],
        [0.4, -0.2, 0.1, -0.51099806603, 0.29961935521, -0.00979594206],
        [0.4, 0.2, -0.1, -0.21849813344, -0.16098019809, 0.16093849962],
    ])
    ref[:, :3] *= angstrom
    check_grid_gradient('co_ccpv5z_pure_hf_g03_fchk', ref, 1e-4)


def check_grid_esp(fn, ref, eps):
    points = ref[:, :3].copy()
    obasis = load_obasis(fn)
    mol = load_mdata(fn)
    dm_full = load_dm(fn)
    esps = obasis.compute_grid_esp_dm(dm_full, mol['coordinates'], mol['pseudo_numbers'], points)
    assert abs(esps - ref[:, 3]).max() < eps


def test_grid_lih_321g_hf_esp_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, 0.906151727538],
        [0.1, 0.0, 0.0, 0.891755005233],
        [0.0, 0.1, 0.0, 0.891755005233],
        [0.0, 0.0, 1.0, 1.422294470114],
        [0.4, 0.2, 0.1, 0.796490099689],
    ])
    ref[:, :3] *= angstrom
    check_grid_esp('li_h_3_21G_hf_g09_fchk', ref, 1e-7)


@attr('slow')
def test_grid_co_ccpv5z_cart_hf_esp_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, 10.69443507172],
        [0.1, 0.0, 0.0, 6.43122889229],
        [0.0, 0.1, 0.0, 6.43406765938],
        [0.0, 0.0, 1.0, 0.27023448629],
        [0.4, 0.2, 0.1, 0.82646540602],
        [-0.4, 0.2, 0.1, 0.93595072191],
        [0.4, -0.2, 0.1, 0.83432301119],
        [0.4, 0.2, -0.1, 0.68524674809],
    ])
    ref[:, :3] *= angstrom
    # cubegen output somehow not reliable?
    check_grid_esp('co_ccpv5z_cart_hf_g03_fchk', ref, 1e-3)


@attr('slow')
def test_grid_co_ccpv5z_pure_hf_esp_some_points():
    ref = np.array([  # from cubegen
        [0.0, 0.0, 0.0, 10.69443507172],
        [0.1, 0.0, 0.0, 6.43122889229],
        [0.0, 0.1, 0.0, 6.43406765938],
        [0.0, 0.0, 1.0, 0.27023448629],
        [0.4, 0.2, 0.1, 0.82646540602],
        [-0.4, 0.2, 0.1, 0.93595072191],
        [0.4, -0.2, 0.1, 0.83432301119],
        [0.4, 0.2, -0.1, 0.68524674809],
    ])
    ref[:, :3] *= angstrom
    check_grid_esp('co_ccpv5z_pure_hf_g03_fchk', ref, 1e-5)


def test_grid_two_index_ne():
    fn = 'li_h_3_21G_hf_g09_fchk'
    obasis = load_obasis(fn)
    mol = load_mdata(fn)
    points, weights = generate_molecular_grid(mol['numbers'], mol['coordinates'], 100000)
    dist0 = np.sqrt(((points - mol['coordinates'][0]) ** 2).sum(axis=1))
    dist1 = np.sqrt(((points - mol['coordinates'][1]) ** 2).sum(axis=1))
    pot = -mol['numbers'][0] / dist0 - mol['numbers'][1] / dist1
    na_ana = obasis.compute_nuclear_attraction(mol['coordinates'], mol['pseudo_numbers'])
    na_grid = obasis.compute_grid_density_fock(points, weights, pot)
    # compare grid-based operator with analytical result
    assert abs(na_grid).max() > 8.0
    np.testing.assert_allclose(na_ana, na_grid, rtol=0.005, atol=0.01)
    # check symmetry
    np.testing.assert_almost_equal(na_grid, na_grid.T)


def test_gob_normalization():
    assert abs(gob_pure_normalization(0.09515, 0) - 0.122100288) < 1e-5
    assert abs(gob_pure_normalization(0.1687144, 1) - 0.154127551) < 1e-5
    assert abs(gob_cart_normalization(0.344, np.array([1, 1, 0])) - 0.440501466) < 1e-8
    assert abs(gob_cart_normalization(0.246, np.array([1, 1, 1])) - 0.242998767) < 1e-8
    assert abs(gob_cart_normalization(0.238, np.array([2, 1, 1])) - 0.127073818) < 1e-8
    assert abs(gob_pure_normalization(0.3, 0) - gob_cart_normalization(0.3, np.array([0, 0, 0]))) < 1e-10
    assert abs(gob_pure_normalization(0.7, 0) - gob_cart_normalization(0.7, np.array([0, 0, 0]))) < 1e-10
    assert abs(gob_pure_normalization(1.9, 0) - gob_cart_normalization(1.9, np.array([0, 0, 0]))) < 1e-10
    assert abs(gob_pure_normalization(0.3, 1) - gob_cart_normalization(0.3, np.array([1, 0, 0]))) < 1e-10
    assert abs(gob_pure_normalization(0.7, 1) - gob_cart_normalization(0.7, np.array([0, 1, 0]))) < 1e-10
    assert abs(gob_pure_normalization(1.9, 1) - gob_cart_normalization(1.9, np.array([0, 0, 1]))) < 1e-10


def test_cart_pure_switch():
    fn = 'water_xyz'
    mol = load_mdata(fn)
    obasis = get_gobasis(mol['coordinates'], mol['numbers'], 'aug-cc-pvdz')
    assert obasis.nbasis == 41
    obasis = get_gobasis(mol['coordinates'], mol['numbers'], 'aug-cc-pvdz', pure=False)
    assert obasis.nbasis == 43


def test_concatenate1():
    fn = 'water_xyz'
    mol = load_mdata(fn)
    obtmp = get_gobasis(mol['coordinates'], mol['numbers'], '3-21G')
    ob = GOBasis.concatenate(obtmp, obtmp)
    assert ob.ncenter == 3 * 2
    assert ob.nbasis == 13 * 2
    a = ob.compute_overlap()
    assert abs(a[:13, :13] - a[:13, 13:]).max() < 1e-15
    assert (a[:13, :13] == a[13:, 13:]).all()
    assert abs(a[:13, :13] - a[13:, :13]).max() < 1e-15


def test_concatenate2():
    fn = 'water_xyz'
    mol = load_mdata(fn)
    obasis1 = get_gobasis(mol['coordinates'], mol['numbers'], '3-21G')
    obasis2 = get_gobasis(mol['coordinates'], mol['numbers'], 'sto-3g')
    obasis = GOBasis.concatenate(obasis1, obasis2)
    assert obasis.ncenter == 3 * 2
    assert obasis.nbasis == obasis1.nbasis + obasis2.nbasis

    a = obasis.compute_overlap()
    a11 = obasis1.compute_overlap()
    a22 = obasis2.compute_overlap()
    N = obasis1.nbasis
    assert (a[:N, :N] == a11).all()
    assert (a[N:, N:] == a22).all()


def test_abstract():
    with assert_raises(NotImplementedError):
        centers = np.zeros((1, 3), float)
        shell_map = np.zeros(2, int)
        nprims = np.array([1, 2])
        shell_types = np.array([0, 1])
        alphas = np.array([1.0, 1.1, 1.2])
        con_coeffs = np.array([0.1, 0.2, 0.3])
        from ..cext import GBasis
        GBasis(centers, shell_map, nprims, shell_types, alphas, con_coeffs)


def test_gobasis_desc_element_map():
    gobd = GOBasisDesc('3-21G', {'H': 'sto-3g', 2: 'cc-pVQZ'})
    coordinates = np.zeros([3, 3])
    numbers = np.array([1, 2, 3])
    obasis = gobd.apply_to(coordinates, numbers)
    assert obasis.centers.shape == (3, 3)
    # H
    assert obasis.shell_map[0] == 0
    assert obasis.nprims[0] == 3
    # He
    assert (obasis.shell_map[1:11] == 1).all()
    assert (obasis.nprims[1:11] == [4, 1, 1, 1, 1, 1, 1, 1, 1, 1]).all()
    # Li
    assert (obasis.shell_map[11:] == 2).all()
    assert (obasis.nprims[11:] == [3, 2, 2, 1, 1]).all()


def test_gobasis_desc_index_map():
    gobd = GOBasisDesc('3-21G', index_map={1: 'sto-3g', 2: 'cc-pVQZ'})
    coordinates = np.zeros([3, 3])
    numbers = np.array([1, 1, 1])
    obasis = gobd.apply_to(coordinates, numbers)
    assert obasis.centers.shape == (3, 3)
    # H
    assert (obasis.shell_map[:2] == 0).all()
    assert (obasis.nprims[:2] == [2, 1]).all()
    # He
    assert (obasis.shell_map[2:3] == 1).all()
    assert (obasis.nprims[2:3] == 3).all()
    # Li
    assert (obasis.shell_map[3:] == 2).all()
    assert (obasis.nprims[3:] == [3, 1, 1, 1, 1, 1, 1, 1, 1, 1]).all()


def test_gobasis_output_args_grid_orbitals_exp():
    fn = 'water_hfs_321g_fchk'
    obasis = load_obasis(fn)
    orb_alpha = load_orbsa_coeffs(fn)
    points = np.random.uniform(-5, 5, (100, 3))
    iorbs = np.array([2, 3])
    orbs1 = np.zeros((100, 2), float)
    obasis.compute_grid_orbitals_exp(orb_alpha, points, iorbs, orbs1)
    orbs2 = obasis.compute_grid_orbitals_exp(orb_alpha, points, iorbs)
    assert (orbs1 == orbs2).all()


def test_gobasis_output_args_grid_density_dm():
    fn = 'water_hfs_321g_fchk'
    obasis = load_obasis(fn)
    points = np.random.uniform(-5, 5, (100, 3))
    rhos1 = np.zeros(100, float)
    dm_full = load_dm(fn)
    obasis.compute_grid_density_dm(dm_full, points, rhos1)
    rhos2 = obasis.compute_grid_density_dm(dm_full, points)
    assert (rhos1 == rhos2).all()


def test_gobasis_output_args_grid_gradient_dm():
    fn = 'water_hfs_321g_fchk'
    obasis = load_obasis(fn)
    points = np.random.uniform(-5, 5, (100, 3))
    gradrhos1 = np.zeros((100, 3), float)
    dm_full = load_dm(fn)
    obasis.compute_grid_gradient_dm(dm_full, points, gradrhos1)
    gradrhos2 = obasis.compute_grid_gradient_dm(dm_full, points)
    assert (gradrhos1 == gradrhos2).all()


def test_gobasis_output_args_grid_hartree_dm():
    fn = 'water_hfs_321g_fchk'
    obasis = load_obasis(fn)
    points = np.random.uniform(-5, 5, (100, 3))
    pots1 = np.zeros(100, float)
    dm_full = load_dm(fn)
    obasis.compute_grid_hartree_dm(dm_full, points, pots1)
    pots2 = obasis.compute_grid_hartree_dm(dm_full, points)
    assert (pots1 == pots2).all()


def test_subset_simple():
    fn = 'water_hfs_321g_fchk'
    obasis = load_obasis(fn)
    # select a basis set for the first hydrogen atom
    sub_obasis, ibasis_list = obasis.get_subset([0, 1])
    assert sub_obasis.ncenter == 1
    assert sub_obasis.nshell == 2
    assert (sub_obasis.centers[0] == obasis.centers[0]).all()
    assert (sub_obasis.shell_map == obasis.shell_map[:2]).all()
    assert (sub_obasis.nprims == obasis.nprims[:2]).all()
    assert (sub_obasis.shell_types == obasis.shell_types[:2]).all()
    assert sub_obasis.nprim_total == 3
    assert (sub_obasis.alphas == obasis.alphas[:3]).all()
    assert (sub_obasis.con_coeffs == obasis.con_coeffs[:3]).all()
    assert (ibasis_list == [0, 1]).all()


def test_subset_simple_reverse():
    fn = 'water_hfs_321g_fchk'
    obasis = load_obasis(fn)
    # select a basis set for the first hydrogen atom
    sub_obasis, ibasis_list = obasis.get_subset([1, 0])
    assert sub_obasis.ncenter == 1
    assert sub_obasis.nshell == 2
    assert (sub_obasis.centers[0] == obasis.centers[0]).all()
    assert (sub_obasis.shell_map == obasis.shell_map[1::-1]).all()
    assert (sub_obasis.nprims == obasis.nprims[1::-1]).all()
    assert (sub_obasis.shell_types == obasis.shell_types[1::-1]).all()
    assert sub_obasis.nprim_total == 3
    assert (sub_obasis.alphas[:1] == obasis.alphas[2:3]).all()
    assert (sub_obasis.alphas[1:] == obasis.alphas[:2]).all()
    assert (sub_obasis.con_coeffs[:1] == obasis.con_coeffs[2:3]).all()
    assert (sub_obasis.con_coeffs[1:] == obasis.con_coeffs[:2]).all()
    assert (ibasis_list == [1, 0]).all()


def test_subset():
    fn = 'water_hfs_321g_fchk'
    obasis = load_obasis(fn)
    # select a basis set for the first hydrogen atom
    sub_obasis, ibasis_list = obasis.get_subset([7, 3, 4, 8])
    assert sub_obasis.ncenter == 2
    assert sub_obasis.nshell == 4
    assert (sub_obasis.centers[0] == obasis.centers[1]).all()
    assert (sub_obasis.centers[1] == obasis.centers[2]).all()
    assert (sub_obasis.shell_map == obasis.shell_map[[7, 3, 4, 8]] - 1).all()
    assert (sub_obasis.nprims == obasis.nprims[[7, 3, 4, 8]]).all()
    assert (sub_obasis.shell_types == obasis.shell_types[[7, 3, 4, 8]]).all()
    assert sub_obasis.nprim_total == 7
    for b0, e0, b1, e1 in (12, 14, 0, 2), (6, 8, 2, 4), (8, 10, 4, 6), (14, 15, 6, 7):
        assert (sub_obasis.alphas[b1:e1] == obasis.alphas[b0:e0]).all()
        assert (sub_obasis.con_coeffs[b1:e1] == obasis.con_coeffs[b0:e0]).all()
    assert (ibasis_list == [11, 3, 4, 5, 6, 12]).all()


def test_basis_atoms():
    fn = 'water_hfs_321g_fchk'
    obasis = load_obasis(fn)
    mol = load_mdata(fn)
    basis_atoms = obasis.get_basis_atoms(mol['coordinates'])
    assert len(basis_atoms) == 3
    icenter = 0
    ibasis_all = []
    for sub_obasis, ibasis_list in basis_atoms:
        assert sub_obasis.ncenter == 1
        assert (sub_obasis.centers[0] == obasis.centers[icenter]).all()
        icenter += 1
        ibasis_all.extend(ibasis_list)
    assert ibasis_all == list(range(obasis.nbasis))


def check_normalization(number, basis):
    """Helper function to test the normalization of contracted basis sets.

    Parameters
    ----------
    number : int
             Element to test. (Keep in mind that not all elements are supported in most
             basis sets.)
    basis : str
            The basis set, e.g. cc-pvdz.
    """
    # Run test on a Helium atom
    coordinates = np.array([[0.0, 0.0, 0.0]])
    numbers = np.array([number])

    # Create a Gaussian basis set
    obasis = get_gobasis(coordinates, numbers, basis)

    # Compute Gaussian integrals
    olp = obasis.compute_overlap()
    np.testing.assert_almost_equal(np.diag(olp), 1.0)


def test_normalization_ccpvdz():
    for number in range(1, 18 + 1):
        check_normalization(number, 'cc-pvdz')
