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


def test_random_rotation():
    for i in xrange(10):
        rotmat = get_random_rotation()
        assert abs(np.dot(rotmat, rotmat.T) - np.identity(3)).max() < 1e-10
        assert abs(np.dot(rotmat.T, rotmat) - np.identity(3)).max() < 1e-10


def test_sphere_int1():
    llgrid = LebedevLaikovSphereGrid(np.zeros(3, float), 1.0, 110)
    x, y, z = llgrid.points.T
    integral = llgrid.integrate(1+x-y+z)
    assert abs(integral - 4*np.pi) < 1e-10


def test_sphere_int2():
    llgrid = LebedevLaikovSphereGrid(np.zeros(3, float), 2.0, 110)
    x, y, z = llgrid.points.T
    integral = llgrid.integrate(1+x-y+z)
    assert abs(integral - 4*np.pi*4) < 1e-10


def test_sphere_int3():
    center = np.random.uniform(-1, 1, 3)
    llgrid = LebedevLaikovSphereGrid(center, 2.0, 110)
    x, y, z = llgrid.points.T
    integral = llgrid.integrate(1+x-y+z)
    check = 4*np.pi*4*(1+center[0]-center[1]+center[2])
    assert abs(integral - check) < 1e-10

def test_sphere_attrs():
    center = np.random.uniform(-1, 1, 3)
    llgrid = LebedevLaikovSphereGrid(center, 2.33, 110)
    assert llgrid.size == 110
    assert llgrid.points.shape == (110, 3)
    assert llgrid.weights.shape == (110,)
    assert llgrid.subgrids is None
    assert (llgrid.center == center).all()
    assert llgrid.radius == 2.33
    assert llgrid.nll == 110
    assert llgrid.random_rotate
