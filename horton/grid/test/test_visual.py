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


def test_line_grid():
    p1 = np.random.normal(0, 1, 3)
    p2 = np.random.normal(0, 1, 3)
    lg = LineGrid(p1, p2, 100)
    assert abs(lg.x[0] + np.linalg.norm(p2-p1)/2) < 1e-10
    assert abs(lg.x[-1] - np.linalg.norm(p2-p1)/2) < 1e-10
    assert abs(lg.points[0] - p1).max() < 1e-10
    assert abs(lg.points[-1] - p2).max() < 1e-10
    lg = LineGrid(p1, p2, 7, extend=1)
    assert abs(lg.x[2] + np.linalg.norm(p2-p1)/2) < 1e-10
    assert abs(lg.x[-3] - np.linalg.norm(p2-p1)/2) < 1e-10
    assert abs(lg.points[2] - p1).max() < 1e-10
    assert abs(lg.points[-3] - p2).max() < 1e-10


def test_rectangle_grid():
    center = np.array([1.0, 0.5, 0.05])
    axis0 = np.array([0.0, 0.0, 0.1])
    axis1 = np.array([0.0, 0.1, 0.0])
    g = RectangleGrid(center, axis0, axis1, -10, 10, -10, 10)
    assert abs(g.weights - 0.01).max() < 1e-10
    assert abs(g.points.mean(axis=0) - center).max() < 1e-10
    assert abs(g.points[0] - np.array([1.0, -0.5, -0.95])).max() < 1e-10
    foo = np.random.normal(0, 1, g.size)
    x, y, z = g.prepare_contour(foo)
    assert x.shape == (21,)
    assert y.shape == (21,)
    assert z.shape == (21, 21)
    assert abs(z.ravel() - foo).max() == 0.0
    assert abs(g.integrate(foo) - foo.sum()*0.01) < 1e-10

    axis0 = np.array([0.0, 0.0, 0.1])
    axis1 = np.array([0.0, 0.1, 0.1])
    g = RectangleGrid(center, axis0, axis1, -10, 10, -10, 10)
    assert abs(g.weights - 0.01).max() < 1e-10
