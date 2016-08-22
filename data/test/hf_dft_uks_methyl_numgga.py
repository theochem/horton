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
    ref_result_kin = 39.32500849927103
    ref_result_energy = -39.762651856264952
    ref_result_grid = 21.728398873370605
    ref_result_exp_alpha = array([-9.91380334, -0.59565247, -0.35562389, -0.35555488, -0.18300771,
        0.06534021,  0.13496298,  0.1349842 ,  0.49468961,  0.53203106,
        0.532239  ,  0.64112016,  0.81382566,  0.81419275,  0.87602631,
        1.58316014,  1.58338864,  1.87853278,  2.03442812,  2.0345952 ])
    ref_result_exp_beta = array([-9.89965882, -0.55986863, -0.34610232, -0.34601783, -0.07913056,
        0.08354667,  0.14220525,  0.14223013,  0.54468578,  0.54489362,
        0.57453112,  0.66511995,  0.81738389,  0.81772358,  0.91285429,
        1.66316173,  1.66320662,  1.98828264,  2.04715243,  2.04734206])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.89584417157022

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
