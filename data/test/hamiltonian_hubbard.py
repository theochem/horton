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


import sys

from numpy import array, allclose
from nose.plugins.attrib import attr

from horton.context import context


@attr('regression_check')
def test_regression():
    ref_result_kin = -7.999999999999998
    ref_result_hartree = 5.9999999999999964
    ref_result_energy = -5.0
    ref_result_ex = -2.9999999999999982
    ref_result_orb_energies = array([ -1.00000000e+00,  -2.39480192e-16,   1.28457889e-16,
         2.00000000e+00,   2.00000000e+00,   3.00000000e+00])

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_ex': 1e-08, 'ref_result_orb_energies': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hamiltonian/hubbard.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
