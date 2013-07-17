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
from nose.tools import assert_raises

from horton import *


def test_interpret_atspec():
    rtf = ExpRTransform(0.1, 1e1, 4)
    rgrid = RadialGrid(rtf, TrapezoidIntegrator1D())

    rgrid0, nlls0 = interpret_atspec(1, 1, (rgrid, 6))
    assert rgrid is rgrid0
    assert (nlls0 == [6,6,6,6]).all()

    rgrid1, nlls1 = interpret_atspec(1, 1, (rgrid, [6,14,26,6]))
    assert rgrid is rgrid1
    assert (nlls1 == [6,14,26,6]).all()

    with assert_raises(ValueError):
        rgrid2, nlls2 = interpret_atspec(1, 1, (rgrid, [1,2,3,4]))
    with assert_raises(ValueError):
        rgrid2, nlls2 = interpret_atspec(1, 1, (rgrid, [6,6,6,6,14]))


def test_atgrid_family_load():
    for af in atgrid_families.itervalues():
        af._load()


def test_atgrid_family_contents1():
    rgrid, nlls = atgrid_families['tv-13.7-3'].get(1, 1)
    assert rgrid.rtransform.to_string() == 'PowerRTransform 7.0879993828935345e-06 16.05937640019924 20'
    assert (nlls == np.array([6, 6, 6, 6, 6, 6, 6, 14, 14, 26, 38, 50, 86, 86, 86, 86, 50, 14, 6, 6])).all()


def test_atgrid_family_contents2():
    rgrid, nlls = atgrid_families['tv-13.7-4'].get(6, 6)
    assert rgrid.rtransform.to_string() == 'PowerRTransform 1.3457673140534789e-07 22.920122695678042 41'
    assert (nlls == np.array([6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 14, 14, 14, 14, 14, 26, 26, 38, 50, 50, 50, 74, 86, 110, 170, 194, 302, 170, 110, 110, 110, 110, 110, 86, 86, 86, 86, 50, 50, 14])).all()


def test_interpret_atspec_family():
    rgrid, nlls = interpret_atspec(1, 1, 'tv-13.7-3')
    assert rgrid.rtransform.to_string() == 'PowerRTransform 7.0879993828935345e-06 16.05937640019924 20'
    assert (nlls == np.array([6, 6, 6, 6, 6, 6, 6, 14, 14, 26, 38, 50, 86, 86, 86, 86, 50, 14, 6, 6])).all()


def test_atomic_grid_basics():
    center = np.random.uniform(-1,1,3)
    rtf = ExpRTransform(0.1, 1e1, 4)
    rgrid = RadialGrid(rtf, TrapezoidIntegrator1D())
    nlls = 6
    for random_rotate in True, False:
        ag0 = AtomicGrid(1, 1, center, (rgrid, 6), random_rotate)
        assert abs(ag0.points.mean(axis=0) - center).max() < 1e-10
        assert (ag0.nlls == [6, 6, 6, 6]).all()
        assert ag0.nsphere == 4
        ag1 = AtomicGrid(1, 1, center, (rgrid, [6, 6, 6, 6]), random_rotate)
        assert abs(ag1.points.mean(axis=0) - center).max() < 1e-10
        assert (ag1.nlls == [6, 6, 6, 6]).all()
        assert ag1.nsphere == 4
        assert abs(ag0.weights - ag1.weights).max() < 1e-10
        assert abs(ag0.av_weights - ag1.av_weights).max() < 1e-10
        assert (abs(ag0.points - ag1.points).max() < 1e-10) ^ random_rotate


def get_hydrogen_1s():
    # density of the 1s orbital
    center = np.random.uniform(-1,1,3)
    rtf = BakerRTransform(2e1, 100)
    rgrid = RadialGrid(rtf, CubicIntegrator1D())
    ag = AtomicGrid(1, 1, center, (rgrid, 110), 100)
    distances = np.sqrt(((center - ag.points)**2).sum(axis=1))
    fn = np.exp(-2*distances)/np.pi
    return ag, fn


def get_hydrogen_1pz():
    # density of the 1pz orbital
    center = np.random.uniform(-1,1,3)
    rtf = BakerRTransform(2e1, 100)
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


def test_spherical_decomposition_hydrogen_1s():
    ag, fn = get_hydrogen_1s()
    sa_fns = ag.get_spherical_average(fn, mtype=2)
    sa_check = np.exp(-2*ag.rgrid.radii)/np.pi
    assert abs(sa_fns[:,0] - sa_check).max() < 1e-10
    assert abs(sa_fns[:,1:]).max() < 1e-10


def test_spherical_decomposition_hydrogen_1pz():
    ag, fn = get_hydrogen_1pz()
    sa_fns = ag.get_spherical_average(fn, mtype=2)
    # s
    sa_check = np.exp(-ag.rgrid.radii)/(32.0*np.pi)*(1.0/3.0)*ag.rgrid.radii**2
    assert abs(ag.rgrid.integrate(sa_check) - 1.0) < 1e-3
    assert abs(sa_fns[:,0] - sa_check).max() < 1e-10
    # p
    assert abs(sa_fns[:,1:4]).max() < 1e-10
    # d
    sa_check = np.exp(-ag.rgrid.radii)/(32.0*np.pi)*(2.0/15.0)*ag.rgrid.radii**4
    assert abs(sa_fns[:,4] - sa_check).max() < 1e-10


def test_atgrid_attrs():
    center = np.array([0.7, 0.2, -0.5], float)
    rtf = ExpRTransform(1e-3, 1e1, 50)
    rgrid = RadialGrid(rtf)
    ag = AtomicGrid(3, 3, center, (rgrid, 26))

    assert ag.size == 50*26
    assert ag.points.shape == (50*26, 3)
    assert ag.weights.shape == (50*26,)
    assert ag.av_weights.shape == (50*26,)
    assert ag.subgrids is None
    assert ag.number == 3
    assert (ag.center == center).all()
    assert ag.rgrid.rtransform == rtf
    assert (ag.nlls == [26]*50).all()
    assert ag.nsphere == 50
    assert ag.random_rotate


def test_random_rotation():
    for i in xrange(10):
        rotmat = get_random_rotation()
        assert abs(np.dot(rotmat, rotmat.T) - np.identity(3)).max() < 1e-10
        assert abs(np.dot(rotmat.T, rotmat) - np.identity(3)).max() < 1e-10
