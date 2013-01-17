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


def test_interpret_atspec():
    rtf = ExpRTransform(0.1, 1e1, 4)
    int1d = TrapezoidIntegrator1D()

    rtf0, int1d0, nlls0 = interpret_atspec(1, (rtf, int1d, 6))
    assert rtf is rtf0
    assert int1d is int1d0
    assert (nlls0 == [6,6,6,6]).all()

    rtf1, int1d1, nlls1 = interpret_atspec(1, (rtf, int1d, [6,14,26,6]))
    assert rtf is rtf1
    assert int1d is int1d1
    assert (nlls1 == [6,14,26,6]).all()

    try:
        rtf2, int1d2, nlls2 = interpret_atspec(1, (rtf, int1d, [1,2,3,4]))
        assert False
    except ValueError:
        pass
    try:
        rtf2, int1d2, nlls2 = interpret_atspec(1, (rtf, int1d, [6,6,6,6,14]))
        assert False
    except ValueError:
        pass


def test_atgrid_family_load():
    for af in atgrid_families.itervalues():
        af._load()


def test_atgrid_family_contents():
    rtf, int1d, nlls = atgrid_families['tv_2012_01_l3'].get(1)
    assert rtf.to_string() == 'ExpRTransform 0.008808017113271033 74.98497608817966 25'
    assert (nlls == np.array([6, 6, 6, 6, 6, 6, 6, 6, 6, 14, 14, 38, 50, 50, 50, 50, 50, 86, 86, 86, 50, 50, 6, 6, 6,])).all()


def test_interpret_atspec_family():
    rtf, int1d, nlls = interpret_atspec(1, 'tv_2012_01_l3')
    assert rtf.to_string() == 'ExpRTransform 0.008808017113271033 74.98497608817966 25'
    assert (nlls == np.array([6, 6, 6, 6, 6, 6, 6, 6, 6, 14, 14, 38, 50, 50, 50, 50, 50, 86, 86, 86, 50, 50, 6, 6, 6,])).all()


def test_atomic_grid_basics():
    center = np.random.uniform(-1,1,3)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(0.1, 1e1, 4)
    nlls = 6
    for random_rotate in True, False:
        ag0 = AtomicGrid(1, center, (rtf, int1d, 6), random_rotate)
        assert abs(ag0.points.mean(axis=0) - center).max() < 1e-10
        assert (ag0.nlls == [6, 6, 6, 6]).all()
        assert ag0.nsphere == 4
        ag1 = AtomicGrid(1, center, (rtf, int1d, [6, 6, 6, 6]), random_rotate)
        assert abs(ag1.points.mean(axis=0) - center).max() < 1e-10
        assert (ag1.nlls == [6, 6, 6, 6]).all()
        assert ag1.nsphere == 4
        assert abs(ag0.weights - ag1.weights).max() < 1e-10
        assert (abs(ag0.points - ag1.points).max() < 1e-10) ^ random_rotate


def test_integrate_hydrogen_1s():
    center = np.random.uniform(-1,1,3)
    int1d = CubicIntegrator1D()
    rtf = BakerRTransform(2e1, 100)
    ag = AtomicGrid(1, center, (rtf, int1d, 110), 100)
    distances = np.sqrt(((center - ag.points)**2).sum(axis=1))
    fn = np.exp(-2*distances)/np.pi
    occupation = ag.integrate(fn)
    assert abs(occupation - 1.0) < 1e-10


def test_atgrid_attrs_1_subgrid():
    center = np.array([0.7, 0.2, -0.5], float)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 50)
    ag = AtomicGrid(2, center, (rtf, int1d, 26), keep_subgrids=1)

    assert ag.size == 50*26
    assert ag.points.shape == (50*26, 3)
    assert ag.weights.shape == (50*26,)
    assert len(ag.subgrids) == 50
    assert ag.number == 2
    assert (ag.center == center).all()
    assert ag.rtransform == rtf
    assert (ag.nlls == [26]*50).all()
    assert ag.nsphere == 50
    assert ag.random_rotate

    radii = rtf.get_radii()
    for j in xrange(50):
        llgrid = ag.subgrids[j]
        assert isinstance(llgrid, LebedevLaikovSphereGrid)
        assert llgrid.size == 26
        assert llgrid.points.shape == (26, 3)
        assert llgrid.weights.shape == (26,)
        assert llgrid.subgrids is None
        assert (llgrid.center == center).all()
        assert llgrid.radius == radii[j]
        assert llgrid.nll == 26
        assert llgrid.random_rotate


def test_atgrid_attrs_0_subgrid():
    center = np.array([0.7, 0.2, -0.5], float)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 50)
    ag = AtomicGrid(3, center, (rtf, int1d, 26), keep_subgrids=0)

    assert ag.size == 50*26
    assert ag.points.shape == (50*26, 3)
    assert ag.weights.shape == (50*26,)
    assert ag.subgrids is None
    assert ag.number == 3
    assert (ag.center == center).all()
    assert ag.rtransform == rtf
    assert (ag.nlls == [26]*50).all()
    assert ag.nsphere == 50
    assert ag.random_rotate
