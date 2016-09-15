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
    ref_result_kin = 39.23271174907517
    ref_result_energy = -39.413743065016263
    ref_result_grid = 22.0721994555369
    ref_result_exp_alpha = array([-9.79899792, -0.59046493, -0.35734025, -0.35723799, -0.17902795,
        0.06451253,  0.12854612,  0.12869712,  0.48639239,  0.51798942,
        0.5182492 ,  0.64737458,  0.81790619,  0.81816594,  0.87705777,
        1.59309677,  1.59333392,  1.89006883,  2.03582521,  2.03683014])
    ref_result_exp_beta = array([-9.78377305, -0.55789987, -0.34289318, -0.34278492, -0.10294034,
        0.07982695,  0.13846989,  0.1386221 ,  0.53416763,  0.53443439,
        0.54845851,  0.66697315,  0.82582109,  0.82606371,  0.90741994,
        1.65792753,  1.65816425,  1.96567071,  2.05634294,  2.05734139])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.79843921229197

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
