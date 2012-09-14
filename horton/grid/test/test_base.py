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


def test_grid_integrate():
    npoint = 10
    grid = BaseGrid(np.random.normal(0, 1, (npoint,3)), np.random.normal(0, 1, npoint))
    pot = np.random.normal(0, 1, npoint)
    dens = np.random.normal(0, 1, npoint)

    # three
    int1 = grid.integrate(pot, dens)
    int2 = (grid.weights*pot*dens).sum()
    assert abs(int1 - int2) < 1e-10

    # two
    int1 = grid.integrate(pot)
    int2 = (grid.weights*pot).sum()
    assert abs(int1 - int2) < 1e-10

    # one
    int1 = grid.integrate()
    int2 = (grid.weights).sum()
    assert abs(int1 - int2) < 1e-10

def test_dot_multi():
    npoint = 10
    pot = np.random.normal(0, 1, npoint)
    dens = np.random.normal(0, 1, npoint)

    # two
    dot1 = np.dot(pot, dens)
    dot2 = dot_multi(pot, dens)
    assert abs(dot1 - dot2) < 1e-10

    # one
    dot1 = pot.sum()
    dot2 = dot_multi(pot)
    assert abs(dot1 - dot2) < 1e-10

    # zero
    dot2 = dot_multi()
    assert dot2 == 0.0
