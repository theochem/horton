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
    ref_result_kin = 39.23193818346811
    ref_result_hartree = 28.075914222094077
    ref_result_energy = -39.412936249624366
    ref_result_grid = -6.002690316053279
    ref_result_exp_alpha = array([-9.79866384, -0.59037636, -0.35713243, -0.35712092, -0.17897514,
        0.06442144,  0.12860056,  0.12861324,  0.48643203,  0.51813855,
        0.51816118,  0.64796082,  0.81865495,  0.81867154,  0.87718107,
        1.59326281,  1.59329876,  1.88965591,  2.0355595 ,  2.03567312])
    ref_result_exp_beta = array([-9.78344212, -0.55780725, -0.34267672, -0.34266459, -0.10290349,
        0.07972093,  0.13851934,  0.13853205,  0.53433239,  0.53435698,
        0.54849195,  0.66758614,  0.82653775,  0.82655137,  0.90749596,
        1.6580803 ,  1.65811823,  1.96525892,  2.05606412,  2.05618171])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.79788328179691

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_beta': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/uks_methyl_lda.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
