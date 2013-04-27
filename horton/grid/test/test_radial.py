# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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


from horton import *

def test_basics1():
    rtf = ExpRTransform(0.1, 1e1, 4)
    int1d = TrapezoidIntegrator1D()
    grid = RadialIntGrid(rtf, int1d)

    assert grid.size == 4
    assert grid.shape == (4,)
    assert grid.rtransform == rtf
    assert grid.int1d == int1d
    assert (grid.weights > 0).all()
    assert (grid.radii == rtf.get_radii()).all()
    assert grid.zeros().shape == (4,)


def test_basics2():
    rtf = BakerRTransform(1e1, 100)
    grid = RadialIntGrid(rtf)

    assert grid.size == 100
    assert grid.shape == (100,)
    assert grid.rtransform == rtf
    assert isinstance(grid.int1d, SimpsonIntegrator1D)
    assert (grid.weights > 0).all()
    assert (grid.radii == rtf.get_radii()).all()
    assert grid.zeros().shape == (100,)


def test_integrate_gauss():
    rtf = ExpRTransform(0.01, 1e1, 100)
    grid = RadialIntGrid(rtf)

    y = np.exp(-0.5*grid.radii**2)
    print grid.integrate(y), (2*np.pi)**1.5
    assert abs(grid.integrate(y) - (2*np.pi)**1.5) < 1e-3
