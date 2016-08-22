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
    ref_result_kin = 39.3491565055777
    ref_result_hartree = 28.101090806446003
    ref_result_energy = -39.829378669752941
    ref_result_ex = -1.230235907807515
    ref_result_grid = -5.206427107498942
    ref_result_exp_alpha = array([-10.20253096,  -0.6640743 ,  -0.40275364,  -0.40273681,
        -0.2248698 ,   0.08886955,   0.16219901,   0.16220686,
         0.53321165,   0.56768587,   0.56770709,   0.68332785,
         0.86346361,   0.86351785,   0.9207139 ,   1.65655084,
         1.6565665 ,   1.95384476,   2.10949156,   2.10952952])
    ref_result_exp_beta = array([-10.18785624,  -0.62231074,  -0.39385382,  -0.39382668,
        -0.04714981,   0.10635468,   0.16869292,   0.16870422,
         0.57658515,   0.57660375,   0.61216788,   0.7033277 ,
         0.86670762,   0.86674928,   0.95059507,   1.72980049,
         1.7299653 ,   2.04441139,   2.12119724,   2.12127052])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.92274790913383

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_ex': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_beta': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/uks_methyl_hybgga.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
