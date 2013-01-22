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


import numpy as np
from horton import *


def test_rectangle_grid():
    center = np.array([1.0, 0.5, 0.05])
    axis0 = np.array([0.0, 0.0, 0.1])
    axis1 = np.array([0.0, 0.1, 0.0])
    g = RectangleGrid(center, axis0, axis1, -10, 10, -10, 10)
    assert abs(g.weights - 0.01).max() < 1e-10
    assert abs(g.points.mean(axis=0) - center).max() < 1e-10
    assert abs(g.points[0] - np.array([1.0, -0.5, -0.95])).max() < 1e-10

    axis0 = np.array([0.0, 0.0, 0.1])
    axis1 = np.array([0.0, 0.1, 0.1])
    g = RectangleGrid(center, axis0, axis1, -10, 10, -10, 10)
    assert abs(g.weights - 0.01).max() < 1e-10

    # TODO test integration weights with a numerical example
