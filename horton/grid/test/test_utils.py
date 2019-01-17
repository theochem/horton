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


from nose.tools import assert_raises
import numpy as np

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_dot_multi():
    cases = [
        (np.arange(10.0), np.arange(10.0, 20.0)),
        (np.random.uniform(0, 1, 5), np.random.uniform(0, 1, 5)),
    ]
    for a, b in cases:
        d1 = np.dot(a, b)
        d2 = dot_multi(a, b)
        assert abs(d1 - d2) < 1e-10

    with assert_raises(AssertionError):
        dot_multi(np.arange(5.0), np.arange(10.0))
