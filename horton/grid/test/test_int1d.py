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
