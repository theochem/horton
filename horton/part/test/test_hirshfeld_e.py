# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2015 The HORTON Development Team
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
#--
#pylint: skip-file


import numpy as np

from horton import *
from horton.part.test.common import get_proatomdb_cp2k


def test_hebasis():
    padb = get_proatomdb_cp2k()
    numbers = np.array([8, 14, 14, 8, 8])
    hebasis = HEBasis(numbers, padb)
    assert hebasis.get_nbasis() == 9 # 3 for every oxygen atom
    assert hebasis.get_atom_nbasis(0) == 3
    assert hebasis.get_atom_nbasis(1) == 0
    assert hebasis.get_atom_nbasis(2) == 0
    assert hebasis.get_atom_nbasis(3) == 3
    assert hebasis.get_atom_nbasis(4) == 3
    for i in 0, 3, 4:
        assert hebasis.get_lower_bound(i, 0) == -1
        assert hebasis.get_lower_bound(i, 1) == -1
        assert hebasis.get_lower_bound(i, 2) == 0
        assert hebasis.get_constant_lico(i) == {0: 1}
    assert (hebasis.get_initial_propars() == 0).all()
    propars = np.array([0.1, 0.5, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    assert hebasis.get_total_lico(0, propars) == {+2: -0.1, +1: -0.4, 0: 1.2, -1: 0.3}
