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
    ref_result_kin = 39.32329980686723
    ref_result_hartree = 28.095118766717036
    ref_result_energy = -39.760924874937807
    ref_result_grid = -6.364574847355393
    ref_result_exp_alpha = array([-9.91322749, -0.59539129, -0.35521695, -0.3552032 , -0.18283053,
        0.06535476,  0.13501027,  0.13502137,  0.49484147,  0.53225531,
        0.53226997,  0.64196942,  0.81493137,  0.81495372,  0.87633336,
        1.58353343,  1.58354862,  1.87807415,  2.03360246,  2.03372245])
    ref_result_exp_beta = array([-9.899089  , -0.55960154, -0.34567608, -0.34566822, -0.07898469,
        0.08352659,  0.14223983,  0.14225618,  0.54492242,  0.5449348 ,
        0.57469671,  0.66601059,  0.81844762,  0.81846533,  0.91309188,
        1.66331277,  1.66344229,  1.98789792,  2.04635398,  2.04642038])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.89455354383031

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_beta': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/uks_methyl_gga.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
