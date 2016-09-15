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
    ref_result_kin = 75.98705562958949
    ref_result_energy = -76.318531129964043
    ref_result_grid = 37.57781551875193
    ref_result_exp_alpha = array([-18.74119815,  -0.9038045 ,  -0.47283089,  -0.29902842,
        -0.22543597,   0.04697382,   0.13005321,   0.75430285,
         0.79457821,   0.84509239,   0.87430993,   1.02013674,
         1.32093985,   1.67019381,   1.67611011,   1.71716065,
         2.22969752,   2.51380346])
    ref_result_nn = 9.1571750364299866
    ref_result_ne = -199.04057731473546

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08}

    test_path = context.get_fn("examples/hf_dft/rks_water_numgga.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
