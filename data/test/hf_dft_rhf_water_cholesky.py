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
    ref_result_kin = 75.97075998377986
    ref_result_x_hf = -8.97383230515447
    ref_result_hartree = 46.89032686623522
    ref_result_energy = -76.025896285286294
    ref_result_exp_alpha = array([-20.54804728,  -1.33138777,  -0.70683717,  -0.55595763,
        -0.4910802 ,   0.18591763,   0.25527901,   0.80543032,
         0.82856233,   1.16143086,   1.20172284,   1.24829854,
         1.46036069,   1.48722528,   1.70052888,   1.87951997,
         1.90390356,   2.46281622,   2.47535859,   3.25634665,
         3.35531301,   3.47264541,   3.9110159 ,   4.1133159 ])
    ref_result_nn = 9.1571750364299866
    ref_result_ne = -199.0703258665769

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_x_hf': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/rhf_water_cholesky.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
