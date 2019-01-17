# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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


def test_stub_weights():
    int1d = StubIntegrator1D()
    assert (int1d.get_weights(4) == np.array([1.0, 1.0, 1.0, 1.0])).all()
    assert len(int1d.get_weights(0)) == 0


def test_trapezoid_weights():
    int1d = TrapezoidIntegrator1D()
    assert (int1d.get_weights(4) == np.array([0.5, 1.0, 1.0, 0.5])).all()


def test_trapezoid_intlin():
    # It should correctly integrate a linear function, duh.
    int1d = TrapezoidIntegrator1D()
    weights = int1d.get_weights(4)
    a, b = np.random.uniform(-1, 1, 2)
    x = np.arange(4)
    y = a*x+b
    int_numer = np.dot(y, weights)
    int_ana = 0.5*(b + a*3+b)*3
    assert abs(int_numer - int_ana) < 1e-10


def test_sine_integration1():
    y = np.sin(np.arange(0,100)/99.0*np.pi)

    int1d = CubicIntegrator1D()
    weights1 = int1d.get_weights(100)
    int_numer1 = np.dot(weights1, y)*np.pi/99.0
    assert abs(int_numer1 - 2.0) < 3e-9

    int2d = TrapezoidIntegrator1D()
    weights2 = int2d.get_weights(100)
    int_numer2 = np.dot(weights2, y)*np.pi/99.0
    assert abs(int_numer2 - 2.0) < 2e-4

    int3d = SimpsonIntegrator1D()
    weights3 = int3d.get_weights(100)
    int_numer3 = np.dot(weights3, y)*np.pi/99.0
    assert abs(int_numer3 - 2.0) < 2e-8

    assert abs(int_numer2 - 2.0) > abs(int_numer3 - 2.0)
    assert abs(int_numer3 - 2.0) > abs(int_numer1 - 2.0)


def test_sine_integration2():
    y = np.sin(np.arange(0,100)/99.0*np.pi/2)

    int1d = CubicIntegrator1D()
    weights1 = int1d.get_weights(100)
    int_numer1 = np.dot(weights1, y)*np.pi/2/99.0
    assert abs(int_numer1 - 1.0) < 1e-7

    int2d = TrapezoidIntegrator1D()
    weights2 = int2d.get_weights(100)
    int_numer2 = np.dot(weights2, y)*np.pi/2/99.0
    assert abs(int_numer2 - 1.0) < 3e-5

    int3d = SimpsonIntegrator1D()
    weights3 = int3d.get_weights(100)
    int_numer3 = np.dot(weights3, y)*np.pi/2/99.0
    assert abs(int_numer3 - 1.0) < 4e-10

    assert abs(int_numer3 - 1.0) < abs(int_numer1 - 1.0)
    assert abs(int_numer1 - 1.0) < abs(int_numer2 - 1.0)
