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

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_mulliken_operators_water_sto3g():
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    mol = IOData.from_file(fn_fchk)
    operators = get_mulliken_operators(mol.obasis, mol.lf)
    for operator in operators:
        assert operator.is_symmetric()
    dm_full = mol.get_dm_full()
    populations = np.array([operator.contract_two('ab,ba', dm_full) for operator in operators])
    charges = mol.numbers - populations
    assert charges[0] < 0 # oxygen atom
    assert abs(charges.sum()) < 1e-3
