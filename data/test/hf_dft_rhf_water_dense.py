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
    ref_result_kin = 75.44720801802484
    ref_result_x_hf = -8.969259906329114
    ref_result_hartree = 46.82886903477147
    ref_result_energy = -75.585812494951796
    ref_result_exp_alpha = array([-20.4241277 ,  -1.32209521,  -0.69192874,  -0.52710998,
        -0.47641493,   0.26331192,   0.36338901,   1.22097216,
         1.27335931,   1.78508596,   1.8583714 ,   2.02905872,   3.10435766])
    ref_result_nn = 9.1571750364299866
    ref_result_ne = -198.04980467784898

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_x_hf': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/rhf_water_dense.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
