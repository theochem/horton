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
    ref_result_kin = 38.93357262027515
    ref_result_hartree = 27.840836401008165
    ref_result_energy = -39.331221904962412
    ref_result_ex = -6.113904009056378
    ref_result_exp_alpha = array([-11.19497791,  -0.92420112,  -0.55513938,  -0.55513937,
        -0.38934657,   0.25358441,   0.3356648 ,   0.33566482,
         0.93322912,   0.98518835,   0.98518849,   1.10244908,
         1.30326226,   1.30326234,   1.67611921])
    ref_result_exp_beta = array([-11.16903157,  -0.81817737,  -0.53903034,  -0.53903034,
         0.16303091,   0.28378927,   0.348972  ,   0.34897202,
         1.00102765,   1.0010278 ,   1.08361697,   1.10609035,
         1.30666579,   1.30666588,   1.72724071])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.07151185985299

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_ex': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_beta': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/uhf_methyl_dense.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
