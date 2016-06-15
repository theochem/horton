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
from nose.plugins.attrib import attr

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import

from horton.test.common import get_random_cell


@attr('slow')
def test_pair_ewald3d_invariance_rcut():
    np.random.seed(0)
    alpha_scale = 4.5
    gcut_scale = 1.5
    delta = np.random.normal(0, 1, 3)
    delta /= np.linalg.norm(delta)
    cell = get_random_cell(1.0, 3)
    results = []
    for rcut in np.arange(10.0, 20.001, 1.0):
        alpha = alpha_scale/rcut
        gcut = gcut_scale*alpha
        results.append(pair_ewald(delta, cell, rcut, alpha, gcut))
    results = np.array(results)
    assert abs(results - results.mean()).max() < 1e-7
