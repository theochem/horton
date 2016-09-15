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
    ref_result_kin = 75.86628050616007
    ref_result_energy = -75.839956410880902
    ref_result_grid = 38.04168508198201
    ref_result_exp_alpha = array([-18.58175488,  -0.89649903,  -0.47396871,  -0.29940619,
        -0.22836067,   0.0451172 ,   0.12846043,   0.76348455,
         0.80251951,   0.83072451,   0.86449305,   1.00433775,
         1.31618315,   1.67084182,   1.67696357,   1.71860491,
         2.23325627,   2.51954412])
    ref_result_nn = 9.1571750364299866
    ref_result_ne = -198.90509703545297

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08}

    test_path = context.get_fn("examples/hf_dft/rks_water_numlda.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
