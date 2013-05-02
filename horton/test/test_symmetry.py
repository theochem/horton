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
from horton.test.common import get_random_cell


def test_symmetry_attrs():
    generators = np.random.uniform(-1, 1, (5, 3, 4))
    fracs = np.random.uniform(0, 1, (4, 3))
    numbers = np.array([1, 6, 6, 1])
    cell = get_random_cell(10.0, 3)
    s = Symmetry('boo', generators, fracs, numbers, cell)
    assert s.labels == ['H0', 'C1', 'C2', 'H3']
    s = Symmetry('boo', generators, fracs, numbers, cell, ['q', 'w', 'e', 'r'])
    assert s.name == 'boo'
    assert s.generators is generators
    assert s.fracs is s.fracs
    assert s.numbers is numbers
    assert s.cell is cell
