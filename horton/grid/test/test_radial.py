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

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_basics1():
    rtf = ExpRTransform(0.1, 1e1, 4)
    int1d = TrapezoidIntegrator1D()
    grid = RadialGrid(rtf, int1d)

    assert grid.size == 4
    assert grid.shape == (4,)
    assert grid.rtransform == rtf
    assert grid.int1d == int1d
    assert (grid.weights > 0).all()
    assert (grid.radii == rtf.get_radii()).all()
    assert grid.zeros().shape == (4,)


def test_basics2():
    rtf = ExpRTransform(1e-3, 1e1, 100)
    grid = RadialGrid(rtf)

    assert grid.size == 100
    assert grid.shape == (100,)
    assert grid.rtransform == rtf
    assert isinstance(grid.int1d, SimpsonIntegrator1D)
    assert (grid.weights > 0).all()
    assert (grid.radii == rtf.get_radii()).all()
    assert grid.zeros().shape == (100,)


def test_integrate_gauss():
    rtf = PowerRTransform(0.0005, 1e1, 100)
    grid = RadialGrid(rtf)
    assert isinstance(grid.int1d, StubIntegrator1D)

    y = np.exp(-0.5*grid.radii**2)
    assert abs(grid.integrate(y) - (2*np.pi)**1.5) < 1e-9
