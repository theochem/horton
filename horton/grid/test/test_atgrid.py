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


def test_get_grid_sizes():
    size, nlls = get_atomic_grid_size(6, 4)
    assert size == 24
    assert (nlls == [6,6,6,6]).all()
    size, nlls = get_atomic_grid_size([6,6,6,6])
    assert size == 24
    assert (nlls == [6,6,6,6]).all()
    size, nlls = get_atomic_grid_size([6,14,26,6])
    assert size == 52
    assert (nlls == [6,14,26,6]).all()
    try:
        get_atomic_grid_size([1,2,3,4])
        assert False
    except ValueError:
        pass
    try:
        get_atomic_grid_size(6)
        assert False
    except ValueError:
        pass


def test_atomic_grid_basics():
    center = np.random.uniform(-1,1,3)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 0.1, 0.3)
    nlls = 6
    for random_rotate in True, False:
        ag0 = AtomicGrid(center, rtf, 6, 4, random_rotate)
        assert abs(ag0.points.mean(axis=0) - center).max() < 1e-10
        assert (ag0.nlls == [6, 6, 6, 6]).all()
        assert ag0.nsphere == 4
        ag1 = AtomicGrid(center, rtf, [6, 6, 6, 6], None, random_rotate)
        assert abs(ag1.points.mean(axis=0) - center).max() < 1e-10
        assert (ag1.nlls == [6, 6, 6, 6]).all()
        assert ag1.nsphere == 4
        assert abs(ag0.weights - ag1.weights).max() < 1e-10
        assert (abs(ag0.points - ag1.points).max() < 1e-10) ^ random_rotate


def test_integrate_hydrogen_1s():
    center = np.random.uniform(-1,1,3)
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 0.1)
    ag = AtomicGrid(center, rtf, 110, 100)
    distances = np.sqrt(((center - ag.points)**2).sum(axis=1))
    fn = np.exp(-2*distances)/np.pi
    occupation = np.dot(fn, ag.weights)
    assert abs(occupation - 1.0) < 1e-8
