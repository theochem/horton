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
from nose.plugins.attrib import attr
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_exceptions():
    with assert_raises(ValueError):
        grid_fn = GB1DMGridDensityFn(-1)

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

    dsq = np.linalg.norm(center - point)**2
    assert abs(work[0] -  scale0*coeff*np.exp(-alpha*dsq)) < 1e-10


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
    dsq = np.linalg.norm(d)**2
    for i in xrange(3):
        assert abs(work[i] -  scales0[i]*coeff*np.exp(-alpha*dsq)*d[i]) < 1e-10


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
    dsq = np.linalg.norm(d)**2
    for i in xrange(3):
        expected = scales0[i]*coeff0*np.exp(-alpha0*dsq)*d[i] + \
                   scales0[i]*coeff1*np.exp(-alpha1*dsq)*d[i]
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
    dsq = np.linalg.norm(d)**2

    expected = scales0[0]*coeff0*np.exp(-alpha0*dsq)*d[0]*d[0] + \
               scales0[0]*coeff1*np.exp(-alpha1*dsq)*d[0]*d[0]
    assert abs(work_cart[0] - expected) < 1e-10 # xx

    expected = scales0[3]*coeff0*np.exp(-alpha0*dsq)*d[1]*d[1] + \
               scales0[3]*coeff1*np.exp(-alpha1*dsq)*d[1]*d[1]
    assert abs(work_cart[3] - expected) < 1e-10 # yy

    expected = scales0[4]*coeff0*np.exp(-alpha0*dsq)*d[1]*d[2] + \
               scales0[4]*coeff1*np.exp(-alpha1*dsq)*d[1]*d[2]
    assert abs(work_cart[4] - expected) < 1e-10 # yz

    grid_fn.cart_to_pure()
    work_pure = grid_fn.get_work(5)
    assert work_pure.shape == (5,)

    from horton.gbasis.test.test_cartpure import tfs
    assert abs(work_pure - np.dot(tfs[2], work_cart)).max() < 1e-10


def test_density_epsilon():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    dm_full = mol.get_dm_full()
    rho1 = mol.obasis.compute_grid_density_dm(dm_full, grid.points)
    for epsilon in 1e-10, 1e-5, 1e-3, 1e-1:
        rho2 = mol.obasis.compute_grid_density_dm(dm_full, grid.points, epsilon=epsilon)
        mask = (rho1 != rho2)
        assert ((rho1[mask] < epsilon) | (abs(rho1[mask]-rho2[mask]) < epsilon)).all()
        assert ((rho2[mask] == 0.0) | (abs(rho1[mask]-rho2[mask]) < epsilon)).all()


def test_density_functional_deriv():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    dm_full = mol.get_dm_full()
    obasis = mol.obasis

    rtf = ExpRTransform(1e-3, 1e1, 10)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 6), random_rotate=False, mode='keep')
    pot = grid.points[:,2].copy()

    def fun(x):
        dm_full._array[:] = x.reshape(obasis.nbasis, -1)
        f = obasis.compute_grid_density_dm(dm_full, grid.points)
        return 0.5*grid.integrate(f, f, pot)

    def fun_deriv(x):
        dm_full._array[:] = x.reshape(obasis.nbasis, -1)
        result = dm_full.new()
        f = obasis.compute_grid_density_dm(dm_full, grid.points)
        obasis.compute_grid_density_fock(grid.points, grid.weights, pot*f, result)
        return result._array.ravel()

    eps = 1e-4
    x = dm_full._array.copy().ravel()
    dxs = []
    for i in xrange(100):
        dxs.append(np.random.uniform(-eps, +eps, x.shape)*x)

    from horton.test.common import check_delta
    check_delta(fun, fun_deriv, x, dxs)


def check_density_gradient(obasis, dm_full, p0, p1):
    points = np.array([p0, p1])
    d = obasis.compute_grid_density_dm(dm_full, points)
    g = obasis.compute_grid_gradient_dm(dm_full, points)
    f1 = d[0] - d[1]
    f2 = np.dot(p0-p1,g[0]+g[1])/2
    assert abs(f1 - f2) < 1e-3*abs(f1)


def test_density_gradient_n2_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    obasis = mol.obasis
    dm_full = mol.get_dm_full()

    points = np.zeros((1,3), float)
    g = obasis.compute_grid_gradient_dm(dm_full, points)
    assert abs(g).max() < 1e-10

    eps = 1e-4
    check_density_gradient(obasis, dm_full, np.array([0.1, 0.3, 0.2]), np.array([0.1+eps, 0.3, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.3, 0.2]), np.array([-0.1+eps, 0.3, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 0.2]), np.array([-0.1+eps, 0.4, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 1.2]), np.array([-0.1+eps, 0.4, 1.2]))


def test_density_gradient_h3_321g():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    obasis = mol.obasis
    dm_full = mol.get_dm_full()

    eps = 1e-4
    check_density_gradient(obasis, dm_full, np.array([0.1, 0.3, 0.2]), np.array([0.1+eps, 0.3, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.3, 0.2]), np.array([-0.1+eps, 0.3, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 0.2]), np.array([-0.1+eps, 0.4, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 1.2]), np.array([-0.1+eps, 0.4, 1.2]))


def test_density_gradient_co_ccpv5z_cart():
    fn_fchk = context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk')
    mol = IOData.from_file(fn_fchk)
    obasis = mol.obasis
    dm_full = mol.get_dm_full()

    eps = 1e-4
    check_density_gradient(obasis, dm_full, np.array([0.1, 0.3, 0.2]), np.array([0.1+eps, 0.3, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.3, 0.2]), np.array([-0.1+eps, 0.3, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 0.2]), np.array([-0.1+eps, 0.4, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 1.2]), np.array([-0.1+eps, 0.4, 1.2]))


def test_density_gradient_co_ccpv5z_pure():
    fn_fchk = context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk')
    mol = IOData.from_file(fn_fchk)
    obasis = mol.obasis
    dm_full = mol.get_dm_full()

    eps = 1e-4
    check_density_gradient(obasis, dm_full, np.array([0.1, 0.3, 0.2]), np.array([0.1+eps, 0.3, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.3, 0.2]), np.array([-0.1+eps, 0.3, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 0.2]), np.array([-0.1+eps, 0.4, 0.2]))
    check_density_gradient(obasis, dm_full, np.array([-0.1, 0.4, 1.2]), np.array([-0.1+eps, 0.4, 1.2]))


def check_dm_gradient(obasis, dm_full, p0, p1):
    grid_fn = GB1DMGridGradientFn(obasis.max_shell_type)

    gradrhos0 = np.zeros((1,3), float)
    obasis._compute_grid1_dm(dm_full, p0, grid_fn, gradrhos0)
    work0 = grid_fn.get_work(grid_fn.max_nbasis)

    gradrhos1 = np.zeros((1,3), float)
    obasis._compute_grid1_dm(dm_full, p1, grid_fn, gradrhos0)
    work1 = grid_fn.get_work(grid_fn.max_nbasis)

    for i in xrange(len(work0)):
        d1 = work0[i,0] - work1[i,0]
        d2 = np.dot(p0-p1, work0[i,1:]+work1[i,1:])/2
        assert abs(d1-d2) < abs(d1)*1e-3


def test_dm_gradient_n2_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    obasis = mol.obasis
    dm_full = mol.get_dm_full()
    eps = 1e-4
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]), np.array([[-0.1+eps, 0.4, 1.2]]))
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]), np.array([[-0.1, 0.4+eps, 1.2]]))
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]), np.array([[-0.1, 0.4, 1.2+eps]]))


def test_dm_gradient_h3_321g():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    obasis = mol.obasis
    dm_full = mol.get_dm_full()
    eps = 1e-4
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]), np.array([[-0.1+eps, 0.4, 1.2]]))
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]), np.array([[-0.1, 0.4+eps, 1.2]]))
    check_dm_gradient(obasis, dm_full, np.array([[-0.1, 0.4, 1.2]]), np.array([[-0.1, 0.4, 1.2+eps]]))


def test_gradient_functional_deriv():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    obasis = mol.obasis
    dm_full = mol.get_dm_full()

    rtf = ExpRTransform(1e-3, 1e1, 10)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 6), random_rotate=False, mode='keep')
    pot = grid.points[:,2].copy()

    def fun(x):
        dm_full._array[:] = x.reshape(obasis.nbasis, -1)
        f = obasis.compute_grid_gradient_dm(dm_full, grid.points)
        tmp = (f*f).sum(axis=1)
        return 0.5*grid.integrate(tmp, pot)

    def fun_deriv(x):
        dm_full._array[:] = x.reshape(obasis.nbasis, -1)
        result = dm_full.new()
        tmp = obasis.compute_grid_gradient_dm(dm_full, grid.points)
        tmp *= pot.reshape(-1,1)
        obasis.compute_grid_gradient_fock(grid.points, grid.weights, tmp, result)
        return result._array.ravel()

    eps = 1e-4
    x = dm_full._array.copy().ravel()
    dxs = []
    for i in xrange(100):
        tmp = np.random.uniform(-eps, +eps, x.shape)*x
        tmp = (tmp+tmp.T)/2
        dxs.append(tmp)

    from horton.test.common import check_delta
    check_delta(fun, fun_deriv, x, dxs)


def check_orbitals(mol):
    points = np.array([
        [0.1, 0.3, 0.2],
        [-0.1, 0.3, 0.2],
        [-0.1, 0.4, 0.2],
        [-0.1, 0.4, 1.2],
    ])

    # just the standard usage (alpha)
    aiorbs = (mol.exp_alpha.occupations > 0).nonzero()[0]
    nalpha = len(aiorbs)
    dm_alpha = mol.exp_alpha.to_dm()
    ad = mol.obasis.compute_grid_density_dm(dm_alpha, points)
    aos = mol.obasis.compute_grid_orbitals_exp(mol.exp_alpha, points, aiorbs)
    ad_check = (aos**2).sum(axis=1)
    assert (abs(ad - ad_check)/abs(ad) < 1e-3).all()

    if hasattr(mol, 'exp_beta'):
        # just the standard usage (beta)
        biorbs = (mol.exp_beta.occupations > 0).nonzero()[0]
        nbeta = len(biorbs)
        dm_beta = mol.exp_beta.to_dm()
        bd = mol.obasis.compute_grid_density_dm(dm_beta, points)
        bos = mol.obasis.compute_grid_orbitals_exp(mol.exp_beta, points, biorbs)
        bd_check = (bos**2).sum(axis=1)
        assert (abs(bd - bd_check)/abs(bd) < 1e-3).all()
    else:
        bos = aos
        bd_check = ad_check

    # compare with full density
    dm_full = mol.get_dm_full()
    fd = mol.obasis.compute_grid_density_dm(dm_full, points)
    fd_check = ad_check + bd_check
    assert (abs(fd - fd_check)/abs(fd) < 1e-3).all()

    # more detailed usage
    assert aos.shape[1] == (mol.exp_alpha.occupations > 0).sum()
    iorbs_alpha = (mol.exp_alpha.occupations > 0).nonzero()[0]
    import random
    iorbs_alpha1 = np.array(random.sample(iorbs_alpha, len(iorbs_alpha)/2))
    iorbs_alpha2 = np.array([i for i in iorbs_alpha if i not in iorbs_alpha1])
    aos1 = mol.obasis.compute_grid_orbitals_exp(mol.exp_alpha, points, iorbs_alpha1)
    aos2 = mol.obasis.compute_grid_orbitals_exp(mol.exp_alpha, points, iorbs_alpha2)
    assert aos1.shape[1] == len(iorbs_alpha1)
    assert aos2.shape[1] == len(iorbs_alpha2)
    ad_check1 = (aos1**2).sum(axis=1)
    ad_check2 = (aos2**2).sum(axis=1)
    assert (abs(ad - ad_check1 - ad_check2)/abs(ad) < 1e-3).all()


def test_orbitals_n2_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    check_orbitals(IOData.from_file(fn_fchk))

def test_orbitals_h3_321g():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    check_orbitals(IOData.from_file(fn_fchk))

def test_orbitals_co_ccpv5z_cart():
    fn_fchk = context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk')
    check_orbitals(IOData.from_file(fn_fchk))

def test_orbitals_co_ccpv5z_pure():
    fn_fchk = context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk')
    check_orbitals(IOData.from_file(fn_fchk))


def check_dm_kinetic(fn, eps=1e-100):
    mol = IOData.from_file(context.get_fn(fn))
    lf = mol.lf
    obasis = mol.obasis
    dm = mol.get_dm_full()
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'fine', random_rotate=False)

    kin = obasis.compute_kinetic(lf)
    ekin1 = kin.contract_two('ab,ab', dm)
    kindens = obasis.compute_grid_kinetic_dm(dm, grid.points)
    ekin2 = grid.integrate(kindens)
    assert abs(ekin1 - ekin2) < eps

    tmp = kin.new()
    obasis.compute_grid_kinetic_fock(grid.points, grid.weights, np.random.uniform(-1, 1, grid.size), tmp)
    assert tmp.is_symmetric()

    kinn = kin.new()
    obasis.compute_grid_kinetic_fock(grid.points, grid.weights, np.ones(grid.size), kinn)
    assert kinn.is_symmetric()
    ekin3 = kinn.contract_two('ab,ab', dm)
    assert abs(ekin1 - ekin3) < eps


def test_dm_kinetic_n2_sto3g():
    check_dm_kinetic('test/n2_hfs_sto3g.fchk', 1e-3)

def test_dm_kinetic_h3_321g():
    check_dm_kinetic('test/h3_pbe_321g.fchk', 5e-5)

@attr('slow')
def test_dm_kinetic_co_ccpv5z_cart():
    check_dm_kinetic('test/co_ccpv5z_cart_hf_g03.fchk', 4e-4)

@attr('slow')
def test_dm_kinetic_co_ccpv5z_pure():
    check_dm_kinetic('test/co_ccpv5z_pure_hf_g03.fchk', 4e-4)


@attr('slow')
def test_kinetic_functional_deriv():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    obasis = mol.obasis
    dm_full = mol.get_dm_full()

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False, mode='keep')
    pot = grid.points[:,2]

    def fun(x):
        dm_full._array[:] = x.reshape(obasis.nbasis, -1)
        f = obasis.compute_grid_kinetic_dm(dm_full, grid.points)
        return 0.5*grid.integrate(f, f, pot)

    def fun_deriv(x):
        dm_full._array[:] = x.reshape(obasis.nbasis, -1)
        result = dm_full.new()
        tmp = obasis.compute_grid_kinetic_dm(dm_full, grid.points)*pot
        obasis.compute_grid_kinetic_fock(grid.points, grid.weights, tmp, result)
        return result._array.ravel()

    eps = 1e-4
    x = dm_full._array.copy().ravel()
    dxs = []
    for i in xrange(100):
        tmp = np.random.uniform(-eps, +eps, x.shape)*x
        dxs.append(tmp)

    from horton.test.common import check_delta
    check_delta(fun, fun_deriv, x, dxs)
