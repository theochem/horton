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



def test_integrate_hydrogen_single_1s():
    numbers = np.array([1], int)
    coordinates = np.array([[0.0, 0.0, -0.5]], float)
    sys = System(coordinates, numbers)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    atspecs = (rtf, 110, 100)

    mg = BeckeMolGrid(sys, atspecs, random_rotate=False)
    dist0 = np.sqrt(((coordinates[0] - mg.points)**2).sum(axis=1))
    fn = np.exp(-2*dist0)/np.pi
    occupation = mg.integrate(fn)
    assert abs(occupation - 1.0) < 1e-3


def test_integrate_hydrogen_pair_1s():
    numbers = np.array([1, 1], int)
    coordinates = np.array([[0.0, 0.0, -0.5], [0.0, 0.0, 0.5]], float)
    sys = System(coordinates, numbers)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    atspecs = (rtf, 110, 100)

    mg = BeckeMolGrid(sys, atspecs, random_rotate=False)
    dist0 = np.sqrt(((coordinates[0] - mg.points)**2).sum(axis=1))
    dist1 = np.sqrt(((coordinates[1] - mg.points)**2).sum(axis=1))
    fn = np.exp(-2*dist0)/np.pi + np.exp(-2*dist1)/np.pi
    occupation = mg.integrate(fn)
    assert abs(occupation - 2.0) < 1e-3


def test_integrate_hydrogen_trimer_1s():
    numbers = np.array([1, 1, 1], int)
    coordinates = np.array([[0.0, 0.0, -0.5], [0.0, 0.0, 0.5], [0.0, 0.5, 0.0]], float)
    sys = System(coordinates, numbers)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    atspecs = (rtf, 110, 100)

    mg = BeckeMolGrid(sys, atspecs, random_rotate=False)
    dist0 = np.sqrt(((coordinates[0] - mg.points)**2).sum(axis=1))
    dist1 = np.sqrt(((coordinates[1] - mg.points)**2).sum(axis=1))
    dist2 = np.sqrt(((coordinates[2] - mg.points)**2).sum(axis=1))
    fn = np.exp(-2*dist0)/np.pi + np.exp(-2*dist1)/np.pi + np.exp(-2*dist2)/np.pi
    occupation = mg.integrate(fn)
    assert abs(occupation - 3.0) < 1e-3


def test_molgrid_attrs_2_subgrid():
    numbers = np.array([6, 8], int)
    coordinates = np.array([[0.0, 0.2, -0.5], [0.1, 0.0, 0.5]], float)
    sys = System(coordinates, numbers)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    atspecs = (rtf, 110, 100)
    mg = BeckeMolGrid(sys, atspecs, keep_subgrids=2)

    assert mg.size == 2*110*100
    assert mg.points.shape == (mg.size, 3)
    assert mg.weights.shape == (mg.size,)
    assert len(mg.subgrids) == 2
    assert mg.system == sys
    assert len(mg.atspecs) == 2
    assert mg.k == 3
    assert mg.random_rotate

    for i in xrange(2):
        atspec = mg.atspecs[i]
        assert atspec[0] == rtf
        assert atspec[1] == [110]*100
        assert atspec[2] == 100*110

        atgrid = mg.subgrids[i]
        assert isinstance(atgrid, AtomicGrid)
        assert atgrid.size == 100*110
        assert atgrid.points.shape == (100*110, 3)
        assert atgrid.weights.shape == (100*110,)
        assert len(atgrid.subgrids) == 100
        assert (atgrid.center == coordinates[i]).all()
        assert atgrid.rtransform == rtf
        assert (atgrid.nlls == [110]*100).all()
        assert atgrid.nsphere == 100
        assert atgrid.random_rotate

        radii = rtf.get_radii(100)
        for j in xrange(100):
            llgrid = atgrid.subgrids[j]
            assert isinstance(llgrid, LebedevLaikovSphereGrid)
            assert llgrid.size == 110
            assert llgrid.points.shape == (110, 3)
            assert llgrid.weights.shape == (110,)
            assert llgrid.subgrids is None
            assert (llgrid.center == coordinates[i]).all()
            assert llgrid.radius == radii[j]
            assert llgrid.nll == 110
            assert llgrid.random_rotate


def test_molgrid_attrs_1_subgrid():
    numbers = np.array([6, 8], int)
    coordinates = np.array([[0.0, 0.2, -0.5], [0.1, 0.0, 0.5]], float)
    sys = System(coordinates, numbers)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    atspecs = (rtf, 110, 100)
    mg = BeckeMolGrid(sys, atspecs, keep_subgrids=1)

    assert mg.size == 2*110*100
    assert mg.points.shape == (mg.size, 3)
    assert mg.weights.shape == (mg.size,)
    assert len(mg.subgrids) == 2
    assert mg.system == sys
    assert len(mg.atspecs) == 2
    assert mg.k == 3
    assert mg.random_rotate

    for i in xrange(2):
        atspec = mg.atspecs[i]
        assert atspec[0] == rtf
        assert atspec[1] == [110]*100
        assert atspec[2] == 100*110

        atgrid = mg.subgrids[i]
        assert isinstance(atgrid, AtomicGrid)
        assert atgrid.size == 100*110
        assert atgrid.points.shape == (100*110, 3)
        assert atgrid.weights.shape == (100*110,)
        assert atgrid.subgrids is None
        assert (atgrid.center == coordinates[i]).all()
        assert atgrid.rtransform == rtf
        assert (atgrid.nlls == [110]*100).all()
        assert atgrid.nsphere == 100
        assert atgrid.random_rotate


def test_molgrid_attrs_0_subgrid():
    numbers = np.array([6, 8], int)
    coordinates = np.array([[0.0, 0.2, -0.5], [0.1, 0.0, 0.5]], float)
    sys = System(coordinates, numbers)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    atspecs = (rtf, 110, 100)
    mg = BeckeMolGrid(sys, atspecs, keep_subgrids=0)

    assert mg.size == 2*110*100
    assert mg.points.shape == (mg.size, 3)
    assert mg.weights.shape == (mg.size,)
    assert mg.subgrids is None
    assert mg.system == sys
    assert len(mg.atspecs) == 2
    assert mg.k == 3
    assert mg.random_rotate
