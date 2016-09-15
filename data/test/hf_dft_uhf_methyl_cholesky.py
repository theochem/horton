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
    ref_result_kin = 39.22375573361731
    ref_result_hartree = 27.9620669215954
    ref_result_energy = -39.554863031594934
    ref_result_ex = -6.150578955382695
    ref_result_exp_alpha = array([-11.26117394,  -0.92269777,  -0.55507866,  -0.55507864,
        -0.38566618,   0.1800017 ,   0.25887874,   0.25887875,
         0.62056367,   0.62056375,   0.66529306,   0.66650271,
         0.86199476,   0.8619948 ,   1.03853021,   1.17555363,
         1.17555367,   1.58300209,   1.66859376,   1.66859393,
         1.70975342,   1.82081328,   2.06943371,   2.06943372,
         2.27776404,   2.27776418,   2.63352911,   2.71827309,   2.71827325])
    ref_result_exp_beta = array([-11.23533385,  -0.82602661,  -0.54105217,  -0.54105216,
         0.14324602,   0.19874758,   0.26771035,   0.26771037,
         0.62291244,   0.62291252,   0.70499609,   0.81625193,
         0.87777073,   0.87777078,   1.07820231,   1.23913672,
         1.23913676,   1.57689179,   1.67276666,   1.67276682,
         1.78656459,   1.82993087,   2.10181673,   2.10181674,
         2.28096012,   2.28096026,   2.63863545,   2.72817412,   2.72817428])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.66989167408859

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_ex': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_beta': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/uhf_methyl_cholesky.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
