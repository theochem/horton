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

from horton import context


@attr('regression_check')
def test_regression():
    ref_result_kin = 39.23261356152271
    ref_result_energy = -39.414134080128676
    ref_result_grid = 22.071676438546728
    ref_result_exp_alpha = array([-9.79909156, -0.59056686, -0.35748404, -0.35733457, -0.17908896,
        0.06441638,  0.12857233,  0.1286138 ,  0.48634891,  0.51789496,
        0.51815698,  0.64746536,  0.81792317,  0.81821749,  0.87695028,
        1.59306753,  1.59339399,  1.89003004,  2.03575936,  2.03707873])
    ref_result_exp_beta = array([-9.78386608, -0.55800322, -0.34303569, -0.34288267, -0.10299764,
        0.07973067,  0.13849873,  0.13853669,  0.53407338,  0.53434975,
        0.54841558,  0.66707389,  0.82584104,  0.82611106,  0.90731201,
        1.65790382,  1.65822529,  1.96563193,  2.05627884,  2.05758675])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.79820902286176

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_beta': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08}

    test_path = context.get_fn("examples/hf_dft/uks_methyl_numlda.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
