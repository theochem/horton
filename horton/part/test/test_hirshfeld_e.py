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
from horton.part.test.common import get_proatomdb_cp2k


def test_hebasis_bound():
    padb = get_proatomdb_cp2k()
    numbers = np.array([8, 14, 14, 8, 8])
    hebasis = HEBasis(numbers, padb)
    for i in 0, 3, 4:
        assert hebasis.get_lower_bound(i, 0) == -1
        assert hebasis.get_lower_bound(i, 1) == -1
        assert hebasis.get_lower_bound(i, 2) == 0
