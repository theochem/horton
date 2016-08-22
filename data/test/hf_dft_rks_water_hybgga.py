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
    ref_result_kin = 76.03393760407292
    ref_result_x_hf = -1.7920659927784908
    ref_result_hartree = 46.89356089593116
    ref_result_energy = -76.406154402221048
    ref_result_grid = -7.568924686343489
    ref_result_exp_alpha = array([-19.12494208,  -0.99561831,  -0.52934337,  -0.35973405,
        -0.28894823,   0.06821496,   0.1532728 ,   0.80077592,
         0.84959969,   0.89304958,   0.92182942,   1.07450833,
         1.37678538,   1.74057685,   1.7462524 ,   1.78615777,
         2.3058015 ,   2.5943212 ])
    ref_result_nn = 9.1571750364299866
    ref_result_ne = -199.12983725953313

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
