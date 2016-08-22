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
    ref_result_kin = 39.23205726235148
    ref_result_hartree = 28.075981237935842
    ref_result_energy = -39.413015749559158
    ref_result_grid = -6.002799780384327
    ref_result_exp_alpha = array([-9.79866525, -0.59038048, -0.35714109, -0.35713059, -0.178971  ,
        0.06442649,  0.12861129,  0.12861694,  0.48643543,  0.51814762,
        0.5181627 ,  0.64794915,  0.81864764,  0.81868303,  0.87718383,
        1.59329324,  1.59330011,  1.8896701 ,  2.03563578,  2.03568129])
    ref_result_exp_beta = array([-9.78344351, -0.55781234, -0.34268607, -0.34267566, -0.10289837,
        0.07972678,  0.13853063,  0.13853633,  0.53434178,  0.53435789,
        0.54849611,  0.66757653,  0.82653056,  0.82656612,  0.90750068,
        1.65811292,  1.65812213,  1.96527178,  2.05614098,  2.05619132])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.79803941212579

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
