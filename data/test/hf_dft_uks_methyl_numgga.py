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
    ref_result_kin = 39.32567734416193
    ref_result_energy = -39.763529945934557
    ref_result_grid = 21.7278763369791
    ref_result_exp_alpha = array([-9.91389785, -0.59579081, -0.35581689, -0.35574836, -0.18304677,
        0.06528646,  0.13499472,  0.13506124,  0.49462852,  0.53203196,
        0.53214151,  0.64129673,  0.81421491,  0.8144156 ,  0.87596567,
        1.58329704,  1.58353004,  1.87865003,  2.03472795,  2.03517142])
    ref_result_exp_beta = array([-9.8997533 , -0.56000854, -0.34629704, -0.34622137, -0.07910827,
        0.08351136,  0.14223621,  0.14230342,  0.54468843,  0.54480147,
        0.57456606,  0.66535111,  0.81775736,  0.81795344,  0.9128101 ,
        1.66320603,  1.66330548,  1.98857455,  2.04744573,  2.04792752])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.89686856973923

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_beta': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08}

    test_path = context.get_fn("examples/hf_dft/uks_methyl_numgga.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
