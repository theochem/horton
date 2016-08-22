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
    ref_result_kin = 108.71316284694448
    ref_result_x_hf = -13.108920740007107
    ref_result_hartree = 74.76419388482375
    ref_result_energy = -108.95408660509986
    ref_result_ne = -302.93789892354584
    ref_result_exp_alpha = array([-15.68653216, -15.68313568,  -1.47097085,  -0.77419995,
        -0.62619794,  -0.60802327,  -0.60802302,   0.17548864,
         0.1754887 ,   0.594534  ,   0.82020386,   0.872305  ,
         0.87230519,   0.99220658,   1.05096918,   1.0509693 ,
         1.14391978,   1.64289937,   1.75749647,   1.75749647,
         1.88172911,   1.88172912,   2.29876807,   2.29876807,
         2.87413368,   2.99611728,   2.99611729,   3.28312431])
    ref_result_nn = 23.615376326684874

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_x_hf': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/rhf_n2_dense.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
