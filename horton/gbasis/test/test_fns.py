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
from nose.plugins.attrib import attr
from nose.tools import assert_raises

from .lightgrid import generate_molecular_grid, integrate
from .common import (load_obasis, load_dm, load_mdata, load_orbsa_coeffs, load_orbsb_coeffs,
                     load_orbsa_occs, load_orbsb_occs, load_orbsa_dms, load_orbsb_dms, check_delta)
from .. import (GOBasis, GOBasisAtom, GOBasisContraction, GB1DMGridGradientFn, GB1DMGridDensityFn,
                get_gobasis)


def check_functional_deriv(fn, comp, dm_method, fock_method):
    """Test consistency of properties on a grid and implementation of Fock build.

    Parameters
    ----------
    fn : str
        The filename of the molecule to run the test on. It should contain a wavefunction
        and it is assumed to be located in the data directory.
    comp : int
        The component of the potential to be tested. (Potentials may have multiple
        components, e.g. in the case of GGA the gradient has different components.)
    dm_method : function
        Computes from a density matrix the density properties of interest, e.g. density.
        This function takes three arguments: obasis, dm_full and points. (See e.g.
        ``GOBasis.compute_grid_density_dm``)
    fock_method : function
        Computes from a potential on a grid, the fock operator. This function takes five
        arguments: obasis, points, weights, potential and fock. (See e.g.
        ``GOBasis.compute_grid_density_fock``)
    """
    obasis = load_obasis(fn)
    dm_full = load_dm(fn)

    mol = load_mdata(fn)
    points, weights = generate_molecular_grid(mol['numbers'], mol['coordinates'])

    def fun(x):
        """Compute the grid properties from a density matrix."""
        dm_full[:] = x.reshape(obasis.nbasis, -1)
        f = dm_method(obasis, dm_full, points)
        f = f.reshape((weights.size, -1))
        return 0.5 * integrate(weights, f[:, comp], f[:, comp])

    def fun_deriv(x):
        """Compute a Fock matrix from potential data on a grid."""
        dm_full[:] = x.reshape(obasis.nbasis, -1)
        tmp = dm_method(obasis, dm_full, points)
        if tmp.ndim > 1:
            tmp[:, :comp] = 0.0
            tmp[:, comp + 1:] = 0.0
        fock = fock_method(obasis, points, weights, tmp)
        return fock.ravel()

    eps = 1e-4
    x = dm_full.copy().ravel()
    dxs = []
    for _irep in range(100):
        tmp = np.random.uniform(-eps, +eps, x.shape) * x
        dxs.append(tmp)

    check_delta(fun, fun_deriv, x, dxs)


def test_exceptions():
    with assert_raises(ValueError):
        GB1DMGridDensityFn(-1)

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])

    with assert_raises(ValueError):
        grid_fn = GB1DMGridDensityFn(2)
        grid_fn.reset(-3, center, point)
    with assert_raises(ValueError):
        grid_fn = GB1DMGridDensityFn(2)
        grid_fn.reset(3, center, point)


def test_grid_fn_s():
    grid_fn = GB1DMGridDensityFn(0)
    assert grid_fn.nwork == 1
    assert grid_fn.max_shell_type == 0
    assert grid_fn.max_nbasis == 1
    assert grid_fn.dim_work == 1
    assert grid_fn.dim_output == 1

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])
    grid_fn.reset(0, center, point)
    assert grid_fn.shell_type0 == 0

    coeff = 0.3
    alpha = 0.5
    scale0 = 0.7
    grid_fn.add(coeff, alpha, np.array([scale0]))
    work = grid_fn.get_work(1)
    assert work.shape == (1,)

    dsq = np.linalg.norm(center - point) ** 2
    assert abs(work[0] - scale0 * coeff * np.exp(-alpha * dsq)) < 1e-10


def test_grid_fn_p():
    grid_fn = GB1DMGridDensityFn(1)
    assert grid_fn.nwork == 3
    assert grid_fn.max_shell_type == 1
    assert grid_fn.max_nbasis == 3
    assert grid_fn.dim_work == 1
    assert grid_fn.dim_output == 1

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])
    grid_fn.reset(1, center, point)
    assert grid_fn.shell_type0 == 1

    coeff = 0.3
    alpha = 0.5
    scales0 = np.array([0.1, 0.2, 0.7])
    grid_fn.add(coeff, alpha, scales0)
    work = grid_fn.get_work(3)
    assert work.shape == (3,)

    d = point - center
    dsq = np.linalg.norm(d) ** 2
    for i in range(3):
        assert abs(work[i] - scales0[i] * coeff * np.exp(-alpha * dsq) * d[i]) < 1e-10


def test_grid_fn_p_contraction():
    grid_fn = GB1DMGridDensityFn(1)
    assert grid_fn.nwork == 3
    assert grid_fn.max_shell_type == 1
    assert grid_fn.max_nbasis == 3
    assert grid_fn.dim_work == 1
    assert grid_fn.dim_output == 1

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])
    grid_fn.reset(1, center, point)
    assert grid_fn.shell_type0 == 1

    scales0 = np.array([0.1, 0.2, 0.7])
    coeff0 = 0.3
    alpha0 = 0.5
    grid_fn.add(coeff0, alpha0, scales0)
    coeff1 = 0.7
    alpha1 = 0.02
    grid_fn.add(coeff1, alpha1, scales0)

    work = grid_fn.get_work(3)
    assert work.shape == (3,)

    d = point - center
    dsq = np.linalg.norm(d) ** 2
    for i in range(3):
        expected = scales0[i] * coeff0 * np.exp(-alpha0 * dsq) * d[i] + \
                   scales0[i] * coeff1 * np.exp(-alpha1 * dsq) * d[i]
        assert abs(work[i] - expected) < 1e-10


def test_grid_fn_d_contraction():
    grid_fn = GB1DMGridDensityFn(3)
    assert grid_fn.nwork == 10
    assert grid_fn.max_shell_type == 3
    assert grid_fn.max_nbasis == 10
    assert grid_fn.dim_work == 1
    assert grid_fn.dim_output == 1

    center = np.array([-0.1, 0.6, -0.3])
    point = np.array([0.5, -0.2, 0.7])
    grid_fn.reset(-2, center, point)
    assert grid_fn.shell_type0 == -2

    scales0 = np.array([0.1, 0.2, 0.7, 0.6, 0.3, 0.8])
    coeff0 = 0.3
    alpha0 = 0.5
    grid_fn.add(coeff0, alpha0, scales0)
    coeff1 = 0.7
    alpha1 = 0.02
    grid_fn.add(coeff1, alpha1, scales0)

    work_cart = grid_fn.get_work(6)
    assert work_cart.shape == (6,)

    d = point - center
    dsq = np.linalg.norm(d) ** 2

    expected = scales0[0] * coeff0 * np.exp(-alpha0 * dsq) * d[0] * d[0] + \
               scales0[0] * coeff1 * np.exp(-alpha1 * dsq) * d[0] * d[0]
    assert abs(work_cart[0] - expected) < 1e-10  # xx

    expected = scales0[3] * coeff0 * np.exp(-alpha0 * dsq) * d[1] * d[1] + \
               scales0[3] * coeff1 * np.exp(-alpha1 * dsq) * d[1] * d[1]
    assert abs(work_cart[3] - expected) < 1e-10  # yy

    expected = scales0[4] * coeff0 * np.exp(-alpha0 * dsq) * d[1] * d[2] + \
               scales0[4] * coeff1 * np.exp(-alpha1 * dsq) * d[1] * d[2]
    assert abs(work_cart[4] - expected) < 1e-10  # yz

    grid_fn.cart_to_pure()
    work_pure = grid_fn.get_work(5)
    assert work_pure.shape == (5,)

    from .test_cartpure import tfs
    assert abs(work_pure - np.dot(tfs[2], work_cart)).max() < 1e-10


def test_density_epsilon():
    fn = 'n2_hfs_sto3g_fchk'
    mol = load_mdata(fn)
    points, weights = generate_molecular_grid(mol['numbers'], mol['coordinates'])
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)
    rho1 = obasis.compute_grid_density_dm(dm_full, points)
    for epsilon in 1e-10, 1e-5, 1e-3, 1e-1:
        rho2 = obasis.compute_grid_density_dm(dm_full, points, epsilon=epsilon)
        mask = (rho1 != rho2)
        assert ((rho1[mask] < epsilon) | (abs(rho1[mask] - rho2[mask]) < epsilon)).all()
        assert ((rho2[mask] == 0.0) | (abs(rho1[mask] - rho2[mask]) < epsilon)).all()


def test_density_functional_deriv():
    check_functional_deriv('n2_hfs_sto3g_fchk', 0, GOBasis.compute_grid_density_dm,
                           GOBasis.compute_grid_density_fock)


def check_density_gradient(obasis, dm_full, point, eps):
    """Finite difference checker for density gradient.

    Parameters
    ----------
    obasis : GOBasis
        The basis set to use in the test.
    dm_full : DenseTwoIndex
        The spin-summed density matrix.
    point : np.ndarray, shape=(3, ), dtype=float
        The reference point for check_delta.
    eps : float
        The magnitude of the displacements.
    """

    def fun(p):
        """Evaluate the density at point p."""
        return obasis.compute_grid_density_dm(dm_full, np.array([p]))[0]

    def fun_deriv(p):
        """Evaluate the density gradient at point p."""
        return obasis.compute_grid_gradient_dm(dm_full, np.array([p]))[0]

    dpoints = np.random.uniform(-eps, eps, (100, 3))
    check_delta(fun, fun_deriv, point, dpoints)


def test_density_gradient_n2_sto3g():
    fn = 'n2_hfs_sto3g_fchk'
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)

    points = np.zeros((1, 3), float)
    g = obasis.compute_grid_gradient_dm(dm_full, points)
    assert abs(g).max() < 1e-10

    eps = 1e-4
    check_density_gradient(obasis, dm_full, np.array([0.1, 0.3, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.3, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 1.2]), eps)


def test_density_gradient_h3_321g():
    fn = 'h3_pbe_321g_fchk'
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)

    eps = 1e-4
    check_density_gradient(obasis, dm_full, np.array([0.1, 0.3, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.3, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 1.2]), eps)


def test_density_gradient_co_ccpv5z_cart():
    fn = 'co_ccpv5z_cart_hf_g03_fchk'
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)

    eps = 1e-4
    check_density_gradient(obasis, dm_full, np.array([0.1, 0.3, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.3, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 1.2]), eps)


def test_density_gradient_co_ccpv5z_pure():
    fn = 'co_ccpv5z_pure_hf_g03_fchk'
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)

    eps = 1e-4
    check_density_gradient(obasis, dm_full, np.array([0.1, 0.3, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.3, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 0.2]), eps)
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 1.2]), eps)


def check_dm_gradient(obasis, dm_full, p0, p1):
    """Check the density gradient with finite differences.

    Parameters
    ----------
    obasis : GOBasis
        The basis set to use in the test.
    dm_full : DenseTwoIndex
        The spin-summed density matrix.
    p0, p1 : np.ndarray, shape=(3, ), dtype=float
        The two points for the finite difference test.
    """
    grid_fn = GB1DMGridGradientFn(obasis.max_shell_type)

    gradrhos0 = np.zeros((1, 3), float)
    obasis._compute_grid1_dm(dm_full, p0, grid_fn, gradrhos0)
    work0 = grid_fn.get_work(grid_fn.max_nbasis)

    gradrhos1 = np.zeros((1, 3), float)
    obasis._compute_grid1_dm(dm_full, p1, grid_fn, gradrhos1)
    work1 = grid_fn.get_work(grid_fn.max_nbasis)

    for i in range(len(work0)):
        d1 = work0[i, 0] - work1[i, 0]
        d2 = np.dot(p0 - p1, work0[i, 1:] + work1[i, 1:]) / 2
        assert abs(d1 - d2) < abs(d1) * 1e-3


def test_dm_gradient_n2_sto3g():
    fn = 'n2_hfs_sto3g_fchk'
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)

    eps = 1e-4
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]),
                      np.array([[-0.1 + eps, 0.4, 1.2]]))
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]),
                      np.array([[-0.1, 0.4 + eps, 1.2]]))
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]),
                      np.array([[-0.1, 0.4, 1.2 + eps]]))


def test_dm_gradient_h3_321g():
    fn = 'h3_hfs_321g_fchk'
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)

    eps = 1e-4
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]),
                      np.array([[-0.1 + eps, 0.4, 1.2]]))
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]),
                      np.array([[-0.1, 0.4 + eps, 1.2]]))
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]),
                      np.array([[-0.1, 0.4, 1.2 + eps]]))


def check_gradient_functional_deriv(fn, comp):
    """Test consistency of density gradient on grid with Fock build."""
    check_functional_deriv(fn, comp, GOBasis.compute_grid_gradient_dm,
                           GOBasis.compute_grid_gradient_fock)


def test_gradient_functional_deriv_0():
    check_gradient_functional_deriv('n2_hfs_sto3g_fchk', 0)


def test_gradient_functional_deriv_1():
    check_gradient_functional_deriv('n2_hfs_sto3g_fchk', 1)


def test_gradient_functional_deriv_2():
    check_gradient_functional_deriv('n2_hfs_sto3g_fchk', 2)


def check_orbitals(fn):
    obasis = load_obasis(fn)

    """Test if sum of squared orbitals reproduces the density, and other things."""
    points = np.array([
        [0.1, 0.3, 0.2],
        [-0.1, 0.3, 0.2],
        [-0.1, 0.4, 0.2],
        [-0.1, 0.4, 1.2],
    ])

    # just the standard usage (alpha)
    orb_alpha = load_orbsa_coeffs(fn)
    orb_alpha_occs = load_orbsa_occs(fn)
    dm_alpha = load_orbsa_dms(fn)
    aiorbs = (orb_alpha_occs > 0).nonzero()[0]
    ad = obasis.compute_grid_density_dm(dm_alpha, points)
    aos = obasis.compute_grid_orbitals_exp(orb_alpha, points, aiorbs)
    ad_check = (aos ** 2).sum(axis=1)
    assert (abs(ad - ad_check) / abs(ad) < 1e-3).all()

    try:
        # just the standard usage (beta)
        orb_beta = load_orbsb_coeffs(fn)
        orb_beta_occs = load_orbsb_occs(fn)
        dm_beta = load_orbsb_dms(fn)
        biorbs = (orb_beta_occs > 0).nonzero()[0]
        bd = obasis.compute_grid_density_dm(dm_beta, points)
        bos = obasis.compute_grid_orbitals_exp(orb_beta, points, biorbs)
        bd_check = (bos ** 2).sum(axis=1)
        assert (abs(bd - bd_check) / abs(bd) < 1e-3).all()
    except IOError:
        bd_check = ad_check

    # compare with full density
    dm_full = load_dm(fn)
    fd = obasis.compute_grid_density_dm(dm_full, points)
    fd_check = ad_check + bd_check
    assert (abs(fd - fd_check) / abs(fd) < 1e-3).all()

    # more detailed usage
    assert aos.shape[1] == (orb_alpha_occs > 0).sum()
    iorbs_alpha = (orb_alpha_occs > 0).nonzero()[0]
    import random
    iorbs_alpha1 = np.array(random.sample(list(iorbs_alpha), int(len(iorbs_alpha) / 2)))
    iorbs_alpha2 = np.array([i for i in iorbs_alpha if i not in iorbs_alpha1])
    aos1 = obasis.compute_grid_orbitals_exp(orb_alpha, points, iorbs_alpha1)
    aos2 = obasis.compute_grid_orbitals_exp(orb_alpha, points, iorbs_alpha2)
    assert aos1.shape[1] == len(iorbs_alpha1)
    assert aos2.shape[1] == len(iorbs_alpha2)
    ad_check1 = (aos1 ** 2).sum(axis=1)
    ad_check2 = (aos2 ** 2).sum(axis=1)
    assert (abs(ad - ad_check1 - ad_check2) / abs(ad) < 1e-3).all()


def test_orbitals_n2_sto3g():
    fn_fchk = 'n2_hfs_sto3g_fchk'
    check_orbitals(fn_fchk)


def test_orbitals_h3_321g():
    fn_fchk = 'h3_pbe_321g_fchk'
    check_orbitals(fn_fchk)


def test_orbitals_co_ccpv5z_cart():
    fn_fchk = 'co_ccpv5z_cart_hf_g03_fchk'
    check_orbitals(fn_fchk)


def test_orbitals_co_ccpv5z_pure():
    fn_fchk = 'co_ccpv5z_pure_hf_g03_fchk'
    check_orbitals(fn_fchk)


def check_dm_kinetic(fn, npoint, eps):
    mol = load_mdata(fn)
    obasis = load_obasis(fn)
    dm = load_dm(fn)
    points, weights = generate_molecular_grid(mol['numbers'], mol['coordinates'], npoint)

    kin = obasis.compute_kinetic()
    ekin1 = np.einsum('ab,ba', kin, dm)
    kindens = obasis.compute_grid_kinetic_dm(dm, points)
    ekin2 = integrate(weights, kindens)
    np.testing.assert_allclose(ekin1, ekin2, atol=eps)

    tmp = obasis.compute_grid_kinetic_fock(
        points, weights, np.random.uniform(-1, 1, weights.size))
    np.testing.assert_almost_equal(tmp, tmp.T)

    kinn = obasis.compute_grid_kinetic_fock(points, weights, np.ones(weights.size))
    np.testing.assert_almost_equal(kinn, kinn.T)
    ekin3 = np.einsum('ab,ba', kinn, dm)
    np.testing.assert_allclose(ekin1, ekin3, atol=eps)


def test_dm_kinetic_n2_sto3g():
    check_dm_kinetic('n2_hfs_sto3g_fchk', 20000, 1e-2)


def test_dm_kinetic_h3_321g():
    check_dm_kinetic('h3_pbe_321g_fchk', 200000, 1e-3)


@attr('slow')
def test_dm_kinetic_co_ccpv5z_cart():
    check_dm_kinetic('co_ccpv5z_cart_hf_g03_fchk', 10000, 0.2)


@attr('slow')
def test_dm_kinetic_co_ccpv5z_pure():
    check_dm_kinetic('co_ccpv5z_pure_hf_g03_fchk', 10000, 0.2)


@attr('slow')
def test_kinetic_functional_deriv():
    check_functional_deriv('n2_hfs_sto3g_fchk', 0, GOBasis.compute_grid_kinetic_dm,
                           GOBasis.compute_grid_kinetic_fock)


def test_concept_gradient():
    # This is not testing a part of Horten directly but it does test a concept
    # used for the evaluation of derivatives of the density.

    # A) Build alphabetically sorted combinations of X, Y and Z
    moms = [[], ['x', 'y', 'z']]
    for l in range(2, 8):
        curmoms = []
        for alpha in 'xyz':
            for oldmom in moms[l - 1]:
                if alpha <= oldmom[0]:
                    curmoms.append(alpha + oldmom)
        moms.append(curmoms)

    # B) Test the rules to get the position of a polynomial of higher order
    #    by adding X, Y or Z
    for l in range(1, 7):
        # rule for x
        for i in range(len(moms[l])):
            s = 'x' + moms[l][i]
            assert s == moms[l + 1][i]
        # rule for y
        for i in range(len(moms[l])):
            s = ''.join(sorted('y' + moms[l][i]))
            nnotx = sum([alpha != 'x' for alpha in moms[l][i]])
            assert s == moms[l + 1][i + 1 + nnotx]
        # rule for z
        for i in range(len(moms[l])):
            s = moms[l][i] + 'z'
            nnotx = sum([alpha != 'x' for alpha in moms[l][i]])
            assert s == moms[l + 1][i + 2 + nnotx]


def check_gradient_systematic(pure):
    """Test density gradient with finite differences on density.

    Parameters
    ----------
    pure : boolean
        If True, pure functions are used in the test, Cartesian otherwise
    """
    # Create fake basis set.
    alpha = 1.5
    bcs = []
    for shell_type in range(2):
        # not properly normalized. So what.
        bcs.append(GOBasisContraction(shell_type, np.array([alpha, alpha / 2]),
                                      np.array([0.5, 0.5])))
    goba = GOBasisAtom(bcs)
    obasis = get_gobasis(np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]), np.array([1, 1]),
                         goba, pure=pure)

    # create fake dm
    dm = np.zeros((obasis.nbasis, obasis.nbasis))

    # Run derivative tests for each DM matrix element.
    eps = 1e-4
    for ibasis0 in range(obasis.nbasis):
        for ibasis1 in range(ibasis0 + 1):
            dm[ibasis0, ibasis1] = 1.2
            dm[ibasis1, ibasis0] = 1.2
            for _irep in range(5):
                point = np.random.normal(0.0, 1.0, 3)
                check_density_gradient(obasis, dm, point, eps)
            dm[ibasis0, ibasis1] = 0.0
            dm[ibasis1, ibasis0] = 0.0


def test_gradient_systematic_cart():
    check_gradient_systematic(False)


def test_gradient_systematic_pure():
    check_gradient_systematic(True)


def check_density_hessian(obasis, dm_full, point, eps):
    """Finite difference checker for density Hessian.

    Parameters
    ----------
    obasis : GOBasis
        The basis set to use in the test.
    dm_full : DenseTwoIndex
        The spin-summed density matrix.
    point : np.ndarray, shape=(3, ), dtype=float
        The reference point for check_delta.
    eps : float
        The magnitude of the displacements.
    """

    def fun(p):
        """Evaluate the density gradient at point p."""
        return obasis.compute_grid_gradient_dm(dm_full, np.array([p]))[0]

    def fun_deriv(p):
        """Evaluate the density Hessian at point p."""
        row = obasis.compute_grid_hessian_dm(dm_full, np.array([p]))[0]
        result = np.zeros((3, 3), float)
        result[0, 0] = row[0]
        result[0, 1] = row[1]
        result[1, 0] = row[1]
        result[0, 2] = row[2]
        result[2, 0] = row[2]
        result[1, 1] = row[3]
        result[1, 2] = row[4]
        result[2, 1] = row[4]
        result[2, 2] = row[5]
        return result

    dpoints = np.random.uniform(-eps, eps, (100, 3))
    check_delta(fun, fun_deriv, point, dpoints)


def check_hessian_systematic(pure):
    """Test density Hessian with finite differences on density gradient.

    Parameters
    ----------
    pure : boolean
        If True, pure functions are used in the test, Cartesian otherwise
    """
    # Create fake basis set.
    alpha = 1.5
    bcs = []
    for shell_type in range(2):
        # not properly normalized. So what.
        bcs.append(GOBasisContraction(shell_type, np.array([alpha, alpha / 2]),
                                      np.array([0.5, 0.5])))
    goba = GOBasisAtom(bcs)
    obasis = get_gobasis(np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]]), np.array([1, 1]),
                         goba, pure=pure)

    # create fake dm
    dm = np.zeros((obasis.nbasis, obasis.nbasis))

    # Run derivative tests for each DM matrix element.
    eps = 1e-4
    for ibasis0 in range(obasis.nbasis):
        for ibasis1 in range(ibasis0 + 1):
            dm[ibasis0, ibasis1] = 1.2
            dm[ibasis1, ibasis0] = 1.2
            for _irep in range(5):
                point = np.random.normal(0.0, 1.0, 3)
                check_density_hessian(obasis, dm, point, eps)
            dm[ibasis0, ibasis1] = 0.0
            dm[ibasis1, ibasis0] = 0.0


def test_hessian_systematic_cart():
    check_hessian_systematic(False)


def test_hessian_systematic_pure():
    check_hessian_systematic(True)


def check_hessian_functional_deriv(fn, comp):
    """Test consistency of Density Hessian on grid with Fock build."""
    check_functional_deriv(fn, comp, GOBasis.compute_grid_hessian_dm,
                           GOBasis.compute_grid_hessian_fock)


def test_hessian_functional_deriv_0():
    check_hessian_functional_deriv('water_sto3g_hf_g03_fchk', 0)


def test_hessian_functional_deriv_1():
    check_hessian_functional_deriv('water_sto3g_hf_g03_fchk', 1)


def test_hessian_functional_deriv_2():
    check_hessian_functional_deriv('water_sto3g_hf_g03_fchk', 2)


def test_hessian_functional_deriv_3():
    check_hessian_functional_deriv('water_sto3g_hf_g03_fchk', 3)


def test_hessian_functional_deriv_4():
    check_hessian_functional_deriv('water_sto3g_hf_g03_fchk', 4)


def test_hessian_functional_deriv_5():
    check_hessian_functional_deriv('water_sto3g_hf_g03_fchk', 5)


def check_gga_evaluation(fn):
    """Test GGA properties on grid. Comparison to related methods.

    Parameters
    ----------
    fn : str
        Filename with wavefunction to use in the test.
    """
    points = np.random.uniform(-5, 5, (1000, 3))
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)

    for _irep in range(5):
        # combined computation of density and gradient
        gga = obasis.compute_grid_gga_dm(dm_full, points)

        # separate computation of density and gradient
        rho = obasis.compute_grid_density_dm(dm_full, points)
        grad = obasis.compute_grid_gradient_dm(dm_full, points)

        assert np.allclose(rho, gga[:, 0], atol=1e-10)
        assert np.allclose(grad, gga[:, 1:4], atol=1e-10)

        # fill the density matrix with random numbers, symmetrize
        dm_full = np.random.normal(0, 1, dm_full.shape)
        dm_full = (dm_full + dm_full.T) / 2


def test_gga_evaluation_co_ccpv5z_cart():
    fn_fchk = 'co_ccpv5z_cart_hf_g03_fchk'
    check_gga_evaluation(fn_fchk)


def test_gga_evaluation_co_ccpv5z_pure():
    fn_fchk = 'co_ccpv5z_pure_hf_g03_fchk'
    check_gga_evaluation(fn_fchk)


def check_gga_fock(fn):
    """Test GGA Fock build. Comparison to related methods.

    Parameters
    ----------
    fn : str
        Filename with wavefunction to use in the test.
    """
    obasis = load_obasis(fn)

    for _irep in range(5):
        # random integration grid
        points = np.random.uniform(-5, 5, (100, 3))
        weights = np.random.uniform(1, 2, 100)
        # combined functional derivatives of toward density and gradient
        pot = np.random.uniform(-1, 1, (100, 4))

        # combined fock matrix build
        fock1 = obasis.compute_grid_gga_fock(points, weights, pot)

        # separate fock matrix build
        fock2 = obasis.compute_grid_density_fock(points, weights, pot[:, 0])
        obasis.compute_grid_gradient_fock(points, weights, pot[:, 1:4], fock2)

        np.testing.assert_almost_equal(fock1, fock2)


def test_gga_fock_co_ccpv5z_cart():
    fn_fchk = 'co_ccpv5z_cart_hf_g03_fchk'
    check_gga_fock(fn_fchk)


def test_gga_fock_co_ccpv5z_pure():
    fn_fchk = 'co_ccpv5z_pure_hf_g03_fchk'
    check_gga_fock(fn_fchk)


def check_gga_functional_deriv(fn, comp):
    """Test consistency of GGA grid properties and Fock build."""
    check_functional_deriv(fn, comp, GOBasis.compute_grid_gga_dm,
                           GOBasis.compute_grid_gga_fock)


def test_gga_functional_deriv_0():
    check_gga_functional_deriv('water_sto3g_hf_g03_fchk', 0)


def test_gga_functional_deriv_1():
    check_gga_functional_deriv('water_sto3g_hf_g03_fchk', 1)


def test_gga_functional_deriv_2():
    check_gga_functional_deriv('water_sto3g_hf_g03_fchk', 2)


def test_gga_functional_deriv_3():
    check_hessian_functional_deriv('water_sto3g_hf_g03_fchk', 3)


def check_mgga_evaluation(fn):
    """Test MGGA properties on grid. Comparison to related methods.

    Parameters
    ----------
    fn : str
        Filename with wavefunction to use in the test.
    """
    # Tests density and gradient by comparing with separate density and gradient
    # evaluation.
    points = np.random.uniform(-5, 5, (1000, 3))
    dm_full = load_dm(fn)
    obasis = load_obasis(fn)

    for _irep in range(5):
        # combined computation of density and gradient
        mgga = obasis.compute_grid_mgga_dm(dm_full, points)

        # separate computation of density and gradient
        rho = obasis.compute_grid_density_dm(dm_full, points)
        grad = obasis.compute_grid_gradient_dm(dm_full, points)
        hess = obasis.compute_grid_hessian_dm(dm_full, points)
        lapl = hess[:, 0] + hess[:, 3] + hess[:, 5]
        tau = obasis.compute_grid_kinetic_dm(dm_full, points)

        assert np.allclose(rho, mgga[:, 0], atol=1e-10)
        assert np.allclose(grad, mgga[:, 1:4], atol=1e-10)
        assert np.allclose(lapl, mgga[:, 4], atol=1e-10)
        assert np.allclose(tau, mgga[:, 5], atol=1e-10)

        # fill the density matrix with random numbers, symmetrize
        dm_full = np.random.normal(0, 1, dm_full.shape)
        dm_full = (dm_full + dm_full.T) / 2


def test_mgga_evaluation_co_ccpv5z_cart():
    fn_fchk = 'co_ccpv5z_cart_hf_g03_fchk'
    check_mgga_evaluation(fn_fchk)


def test_mgga_evaluation_co_ccpv5z_pure():
    fn_fchk = 'co_ccpv5z_pure_hf_g03_fchk'
    check_mgga_evaluation(fn_fchk)


def check_mgga_fock(fn):
    """Test MGGA Fock build. Comparison to related methods.

    Parameters
    ----------
    fn : str
        Filename with wavefunction to use in the test.
    """
    obasis = load_obasis(fn)

    for _irep in range(5):
        # random integration grid
        points = np.random.uniform(-5, 5, (100, 3))
        weights = np.random.uniform(1, 2, 100)
        # combined functional derivatives of toward density and gradient
        pot = np.random.uniform(-1, 1, (100, 6))

        # combined fock matrix build
        fock1 = obasis.compute_grid_mgga_fock(points, weights, pot)

        # separate fock matrix build
        fock2 = obasis.compute_grid_density_fock(points, weights, pot[:, 0])
        obasis.compute_grid_gradient_fock(points, weights, pot[:, 1:4], fock2)
        hessian_pot = np.zeros((100, 6), float)
        hessian_pot[:, 0] = pot[:, 4]
        hessian_pot[:, 3] = pot[:, 4]
        hessian_pot[:, 5] = pot[:, 4]
        obasis.compute_grid_hessian_fock(points, weights, hessian_pot, fock2)
        obasis.compute_grid_kinetic_fock(points, weights, pot[:, 5], fock2)

        np.testing.assert_almost_equal(fock1, fock2)


def test_mgga_fock_co_ccpv5z_cart():
    fn_fchk = 'co_ccpv5z_cart_hf_g03_fchk'
    check_mgga_fock(fn_fchk)


def test_mgga_fock_co_ccpv5z_pure():
    fn_fchk = 'co_ccpv5z_pure_hf_g03_fchk'
    check_mgga_fock(fn_fchk)


def check_mgga_functional_deriv(fn, comp):
    """Test consistency of MGGA grid properties and Fock build."""
    check_functional_deriv(fn, comp, GOBasis.compute_grid_mgga_dm,
                           GOBasis.compute_grid_mgga_fock)


def test_mgga_functional_deriv_0():
    check_mgga_functional_deriv('water_sto3g_hf_g03_fchk', 0)


def test_mgga_functional_deriv_1():
    check_mgga_functional_deriv('water_sto3g_hf_g03_fchk', 1)


def test_mgga_functional_deriv_2():
    check_mgga_functional_deriv('water_sto3g_hf_g03_fchk', 2)


def test_mgga_functional_deriv_3():
    check_mgga_functional_deriv('water_sto3g_hf_g03_fchk', 3)


def test_mgga_functional_deriv_4():
    check_mgga_functional_deriv('water_sto3g_hf_g03_fchk', 4)


def test_mgga_functional_deriv_5():
    check_mgga_functional_deriv('water_sto3g_hf_g03_fchk', 5)
