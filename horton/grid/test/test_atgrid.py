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


import os, shutil, numpy as np, h5py as h5
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import

from horton.test.common import tmpdir


def test_normalize_nlls():
    from horton.grid.atgrid import _normalize_nlls
    assert (_normalize_nlls(6, 10) == np.array([6]*10)).all()
    assert (_normalize_nlls([6], 10) == np.array([6]*10)).all()
    assert (_normalize_nlls([6, 6, 6], 3) == np.array([6]*3)).all()
    with assert_raises(ValueError):
        _normalize_nlls([6, 6, 6], 4)


def test_agspec_tuple1():
    rtf = ExpRTransform(0.1, 1e1, 4)
    rgrid = RadialGrid(rtf, StubIntegrator1D())

    agspec = AtomicGridSpec((rgrid, 6))
    rgrid0, nlls0 = agspec.get(1, 1)
    assert rgrid is rgrid0
    assert (nlls0 == [6,6,6,6]).all()
    rgrid0, nlls0 = agspec.get(2, 2)
    assert rgrid is rgrid0
    assert (nlls0 == [6,6,6,6]).all()


def test_agspec_tuple2():
    rtf = ExpRTransform(0.1, 1e1, 4)
    rgrid = RadialGrid(rtf, StubIntegrator1D())

    agspec = AtomicGridSpec((rgrid, [6,14,26,6]))
    rgrid0, nlls0 = agspec.get(1, 1)
    assert rgrid is rgrid0
    assert (nlls0 == [6,14,26,6]).all()

    with assert_raises(ValueError):
        agspec = AtomicGridSpec((rgrid, [6,14,26,6,6]))


def test_agspec_list():
    agspec = AtomicGridSpec([
        (1, 1, RadialGrid(ExpRTransform(0.1, 1e1, 4), StubIntegrator1D()), [6,14,26,6]),
        (2, 2, RadialGrid(ExpRTransform(0.2, 1e1, 4), StubIntegrator1D()), [6,14,26,14]),
        (10, 8, RadialGrid(ExpRTransform(0.3, 1e1, 4), StubIntegrator1D()), [6,14,26,26]),
        (10, 10, RadialGrid(ExpRTransform(0.4, 1e1, 4), StubIntegrator1D()), [6,14,26,38]),
    ])
    rgrid, nlls = agspec.get(1, 1)
    assert rgrid.rtransform.rmin == 0.1
    assert (nlls == [6,14,26,6]).all()
    rgrid, nlls = agspec.get(2, 2)
    assert rgrid.rtransform.rmin == 0.2
    assert (nlls == [6,14,26,14]).all()
    rgrid, nlls = agspec.get(10, 8)
    assert rgrid.rtransform.rmin == 0.3
    assert (nlls == [6,14,26,26]).all()
    rgrid, nlls = agspec.get(10, 10)
    assert rgrid.rtransform.rmin == 0.4
    assert (nlls == [6,14,26,38]).all()
    rgrid, nlls = agspec.get(10, 6)
    assert rgrid.rtransform.rmin == 0.3
    assert (nlls == [6,14,26,26]).all()
    rgrid, nlls = agspec.get(2, 1)
    assert rgrid.rtransform.rmin == 0.2
    assert (nlls == [6,14,26,14]).all()


def test_agspec_string():
    agspec = AtomicGridSpec('power:0.001:10.0:20:26')
    for number in 1, 4, 10:
        rgrid, nlls = agspec.get(number, number)
        assert isinstance(rgrid.rtransform, PowerRTransform)
        assert rgrid.rtransform.rmin == 0.001*angstrom
        assert rgrid.rtransform.rmax == 10.0*angstrom
        assert rgrid.size == 20
        assert (nlls == [26]*20).all()


def test_agspec_wrong_string():
    with assert_raises(ValueError):
        AtomicGridSpec('power:0.001:10.0:20')


def test_agspec_local_file():
    with tmpdir('horton.scripts.test.test_espfit.test_scripts_symmetry') as dn:
        fn_dest = os.path.join(dn, 'mygrid.txt')
        shutil.copy(context.get_fn('grids/tv-13.7-4.txt'), fn_dest)
        agspec = AtomicGridSpec(fn_dest)
        rgrid, nlls = agspec.get(1, 1)
        assert rgrid.rtransform.to_string() == 'PowerRTransform 3.69705074304963e-06 19.279558946793685 24'
        assert (nlls == np.array([6, 6, 6, 6, 6, 6, 6, 6, 14, 14, 26, 38, 50, 86, 110, 110, 110, 110, 86, 50, 50, 14, 6, 6])).all()


def test_agspec_load_simple_names():
    nrads = [20, 24, 34, 41, 49, 59]
    for name in 'coarse', 'medium', 'fine', 'veryfine', 'ultrafine', 'insane':
        agspec = AtomicGridSpec(name)
        assert len(agspec.get(1, 1)[1]) == nrads.pop(0)


def test_agspec_load_names():
    nrads = [20, 24, 34, 41, 49, 59]
    for name in 'tv-13.7-3', 'tv-13.7-4', 'tv-13.7-5', 'tv-13.7-6', 'tv-13.7-7', 'tv-13.7-8':
        agspec = AtomicGridSpec(name)
        assert len(agspec.get(1, 1)[1]) == nrads.pop(0)


def test_agspec_coarse_contents():
    rgrid, nlls = AtomicGridSpec('coarse').get(1, 1)
    assert rgrid.rtransform.to_string() == 'PowerRTransform 7.0879993828935345e-06 16.05937640019924 20'
    assert (nlls == np.array([6, 6, 6, 6, 6, 6, 6, 14, 14, 26, 38, 50, 86, 86, 86, 86, 50, 14, 6, 6])).all()


def test_atgrid_medium_contents():
    for agspec in AtomicGridSpec('medium'), AtomicGridSpec():
        rgrid, nlls = agspec.get(1, 1)
        assert rgrid.rtransform.to_string() == 'PowerRTransform 3.69705074304963e-06 19.279558946793685 24'
        assert (nlls == np.array([6, 6, 6, 6, 6, 6, 6, 6, 14, 14, 26, 38, 50, 86, 110, 110, 110, 110, 86, 50, 50, 14, 6, 6])).all()


def test_agspec_get_size():
    agspec = AtomicGridSpec()
    assert agspec.get_size(1, 1) == 928
    assert agspec.get_size(8, 8) == 3754


def test_atomic_grid_basics1():
    center = np.random.uniform(-1,1,3)
    rtf = ExpRTransform(0.1, 1e1, 4)
    rgrid = RadialGrid(rtf, StubIntegrator1D())
    for random_rotate in True, False:
        ag0 = AtomicGrid(1, 1, center, (rgrid, 6), random_rotate)
        assert abs(ag0.points.mean(axis=0) - center).max() < 1e-10
        assert (ag0.nlls == 6).all()
        assert ag0.nsphere == 4
        assert (ag0.lmaxs == 3).all()
        ag1 = AtomicGrid(1, 1, center, (rgrid, [6, 6, 6, 6]), random_rotate)
        assert abs(ag1.points.mean(axis=0) - center).max() < 1e-10
        assert (ag0.nlls == 6).all()
        assert ag0.nsphere == 4
        assert (ag0.lmaxs == 3).all()
        assert abs(ag0.weights - ag1.weights).max() < 1e-10
        assert (abs(ag0.points - ag1.points).max() < 1e-10) ^ random_rotate


def test_atomic_grid_basics2():
    center = np.random.uniform(-1,1,3)
    rtf = ExpRTransform(0.1, 1e1, 3)
    rgrid = RadialGrid(rtf, StubIntegrator1D())
    ag2 = AtomicGrid(1, 1, center, (rgrid, [6, 14, 26]))
    assert abs(ag2.points.mean(axis=0) - center).max() < 1e-10
    assert (ag2.nlls == [6, 14, 26]).all()
    assert ag2.nsphere == 3
    assert (ag2.lmaxs == [3, 5, 7]).all()


def get_hydrogen_1s():
    # density of the 1s orbital
    center = np.random.uniform(-1,1,3)
    rtf = PowerRTransform(1e-3, 2e1, 100)
    rgrid = RadialGrid(rtf, CubicIntegrator1D())
    ag = AtomicGrid(1, 1, center, (rgrid, 110), 100)
    distances = np.sqrt(((center - ag.points)**2).sum(axis=1))
    fn = np.exp(-2*distances)/np.pi
    return ag, fn


def get_hydrogen_1pz():
    # density of the 1pz orbital
    center = np.random.uniform(-1,1,3)
    rtf = PowerRTransform(1e-3, 2e1, 100)
    rgrid = RadialGrid(rtf, SimpsonIntegrator1D())
    ag = AtomicGrid(1, 1, center, (rgrid, 110), 100)
    z = ag.points[:,2] - center[2]
    distances = np.sqrt(((center - ag.points)**2).sum(axis=1))
    fn = np.exp(-distances)/(32.0*np.pi)*z**2
    return ag, fn


def test_integrate_hydrogen_1s():
    ag, fn = get_hydrogen_1s()
    occupation = ag.integrate(fn)
    assert abs(occupation - 1.0) < 1e-10


def test_spherical_average_hydrogen_1s():
    ag, fn = get_hydrogen_1s()
    x = ag.points[:,0] - ag.center[0]
    y = ag.points[:,1] - ag.center[1]
    z = ag.points[:,2] - ag.center[2]
    sa_check = np.exp(-2*ag.rgrid.radii)/np.pi
    for cx, cy, cz, cxxx in (0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 1, 0, 0), (0, 1, 0, 1):
        sa_fn = ag.get_spherical_average(fn + cx*x + cy*y + cz*z + cxxx*x*x*x)
        assert abs(sa_fn - sa_check).max() < 1e-10


def test_spherical_average_grads1():
    center = np.random.uniform(-1,1,3)
    rtf = PowerRTransform(1e-3, 2e1, 100)
    rgrid = RadialGrid(rtf, CubicIntegrator1D())
    ag = AtomicGrid(1, 1, center, (rgrid, 110), 100)

    delta = ag.points - ag.center
    d = np.sqrt(delta[:,0]**2 + delta[:,1]**2 + delta[:,2]**2)

    fy1 = np.exp(-d**2)
    fg1 = -2*delta*(np.exp(-d**2)).reshape(-1,1)

    r = ag.rgrid.rtransform.get_radii()
    say, sad = ag.get_spherical_average(fy1, grads=[fg1])
    say_check = np.exp(-r**2)
    sad_check = -2*r*np.exp(-r**2)
    assert abs(say - say_check).max() < 1e-10
    assert abs(sad - sad_check).max() < 1e-10


def test_spherical_average_grads2():
    center = np.random.uniform(-1,1,3)
    rtf = PowerRTransform(1e-3, 2e1, 100)
    rgrid = RadialGrid(rtf, CubicIntegrator1D())
    ag = AtomicGrid(1, 1, center, (rgrid, 110), 100)

    delta = ag.points - ag.center
    d = np.sqrt(delta[:,0]**2 + delta[:,1]**2 + delta[:,2]**2)

    fy1 = np.exp(-d**2)
    fg1 = -2*delta*(np.exp(-d**2)).reshape(-1,1)
    fy2 = np.exp(-d)
    fg2 = -delta*(np.exp(-d)/d).reshape(-1,1)

    r = ag.rgrid.rtransform.get_radii()
    say, sad = ag.get_spherical_average(fy1, fy2, grads=[fg1, fg2])
    say_check = np.exp(-r**2 - r)
    sad_check = (-2*r - 1)*np.exp(-r**2 - r)
    assert abs(say - say_check).max() < 1e-10
    assert abs(sad - sad_check).max() < 1e-10


def test_spherical_decomposition_hydrogen_1s():
    ag, fn = get_hydrogen_1s()
    sa_fns = ag.get_spherical_decomposition(fn, lmax=4)
    sa_check = np.exp(-2*ag.rgrid.radii)/np.pi*np.sqrt(4*np.pi)
    assert abs(sa_fns[0].y - sa_check).max() < 1e-10
    for sa_fn in sa_fns[1:]:
        assert abs(sa_fn.y).max() < 1e-10


def test_spherical_decomposition_hydrogen_1pz():
    ag, fn = get_hydrogen_1pz()
    sa_fns = ag.get_spherical_decomposition(fn, lmax=4)
    # s
    sa_check = np.exp(-ag.rgrid.radii)/(32.0*np.pi)*(1.0/3.0)*ag.rgrid.radii**2*np.sqrt(4*np.pi)
    assert abs(sa_fns[0].y - sa_check).max() < 1e-10
    # p
    for sa_fn in sa_fns[1:4]:
        assert abs(sa_fn.y).max() < 1e-10
    # d
    sa_check = np.exp(-ag.rgrid.radii)/(32.0*np.pi)*(2.0/15.0)*ag.rgrid.radii**2*np.sqrt(5*4*np.pi)
    assert abs(sa_fns[4].y - sa_check).max() < 1e-10


def test_spherical_decomposition_hydrogen_1pz_conventions():
    ag, fn = get_hydrogen_1pz()
    sa_fns = ag.get_spherical_decomposition(fn, lmax=4)
    multipoles = ag.integrate(fn, center=ag.center, mtype=2, lmax=4)
    # s
    assert abs(ag.rgrid.integrate(sa_fns[0].y) - np.sqrt(4*np.pi)) < 1e-3
    assert abs(multipoles[0] - 1.0) < 1e-3
    # d
    r = ag.rgrid.rtransform.get_radii()
    qzz = ag.rgrid.integrate(sa_fns[4].y, r, r)/np.sqrt(5*4*np.pi)
    print multipoles[4], qzz
    assert abs(multipoles[4] - qzz) < 1e-10


def test_atgrid_attrs():
    center = np.array([0.7, 0.2, -0.5], float)
    rtf = ExpRTransform(1e-3, 1e1, 50)
    rgrid = RadialGrid(rtf)
    ag = AtomicGrid(3, 3, center, (rgrid, 26))

    assert ag.size == 50*26
    assert ag.points.shape == (50*26, 3)
    assert ag.weights.shape == (50*26,)
    assert ag.subgrids is None
    assert ag.number == 3
    assert (ag.center == center).all()
    assert ag.rgrid.rtransform == rtf
    assert (ag.nlls == [26]*50).all()
    assert ag.nsphere == 50
    assert ag.random_rotate


def test_random_rotation():
    for _ in xrange(10):
        rotmat = get_random_rotation()
        assert abs(np.dot(rotmat, rotmat.T) - np.identity(3)).max() < 1e-10
        assert abs(np.dot(rotmat.T, rotmat) - np.identity(3)).max() < 1e-10


def test_agspec_hdf5_coarse():
    agspec1 = AtomicGridSpec('coarse')
    with h5.File('horton.grid.test.test_atgrid.test_agspec_hdf5_coarse', driver='core', backing_store=False) as f:
        agspec1.to_hdf5(f)
        agspec2 = AtomicGridSpec.from_hdf5(f)
    assert sorted(agspec1.members.keys()) == sorted(agspec2.members.keys())
    for number, cases1 in agspec1.members.iteritems():
        cases2 = agspec2.members[number]
        assert len(cases1) == len(cases2)
        for case1, case2 in zip(cases1, cases2):
            assert case1[0] == case2[0]
            assert case1[1].rtransform.to_string() == case2[1].rtransform.to_string()
            assert (case1[2] == case2[2]).all()
