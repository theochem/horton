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
    ref_result_kin = 76.03389826795868
    ref_result_x_hf = -1.7920648352546238
    ref_result_hartree = 46.89351081194541
    ref_result_energy = -76.406157158872574
    ref_result_grid = -7.568922226522652
    ref_result_exp_alpha = array([-19.12495235,  -0.99562081,  -0.52934442,  -0.3597406 ,
        -0.28895222,   0.0681964 ,   0.15329257,   0.80078421,
         0.84959223,   0.89305046,   0.92182136,   1.07452226,
         1.37676279,   1.74056963,   1.74623791,   1.78613112,
         2.30577158,   2.59429795])
    ref_result_nn = 9.1571750364299866
    ref_result_ne = -199.12975421342938

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_x_hf': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/rks_water_hybgga.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
