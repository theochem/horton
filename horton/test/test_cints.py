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


def test_fact():
    assert fact(-20) == 1
    assert fact(0) == 1
    assert fact(1) == 1
    assert fact(2) == 2
    assert fact(3) == 6
    assert fact(4) == 24
    assert fact(5) == 120


def test_fact2():
    assert fact2(-20) == 1
    assert fact2(0) == 1
    assert fact2(1) == 1
    assert fact2(2) == 2
    assert fact2(3) == 3
    assert fact2(4) == 8
    assert fact2(5) == 15


def test_overlap():
    for nx in xrange(3):
        for ny in xrange(3):
            for nz in xrange(3):
                for alpha in np.arange(0.5, 2.51, 0.5):
                    x, y, z = np.random.uniform(-1, 1, 3)
                    olp = cints_overlap(alpha, nx, ny, nz, x, y, z,
                                        alpha, nx, ny, nz, x, y, z)
                    check = olp*gob_normalization(alpha, nx, ny, nz)**2
                    assert abs(check - 1.0) < 1e-10
