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
    ref_result_kin = 39.32330866679178
    ref_result_hartree = 28.095142335211435
    ref_result_energy = -39.760960634478579
    ref_result_grid = -6.364617950299294
    ref_result_exp_alpha = array([-9.91322565, -0.5953944 , -0.35521501, -0.35521421, -0.18283186,
        0.06535657,  0.13501949,  0.13502512,  0.49484186,  0.53224793,
        0.53226908,  0.64196645,  0.81492654,  0.81497286,  0.87633934,
        1.5835337 ,  1.58354158,  1.87809065,  2.03363447,  2.03369493])
    ref_result_exp_beta = array([-9.89908817, -0.55959816, -0.34568878, -0.34567415, -0.07893459,
        0.08353716,  0.14223953,  0.14225655,  0.54491004,  0.54493374,
        0.5747816 ,  0.66601651,  0.81842748,  0.81848177,  0.91311294,
        1.66319503,  1.66341063,  1.9881187 ,  2.04635541,  2.0464332 ])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.89457862884615

    thresholds = {'ref_result_kin': 1e-08, 'ref_result_energy': 1e-08, 'ref_result_grid': 1e-08, 'ref_result_ne': 1e-08, 'ref_result_exp_beta': 1e-08, 'ref_result_exp_alpha': 1e-08, 'ref_result_nn': 1e-08, 'ref_result_hartree': 1e-08}

    test_path = context.get_fn("examples/hf_dft/uks_methyl_gga.py")

    l = {}
    m = locals()
    with open(test_path) as fh:
        exec fh in l

    for k,v in thresholds.items():
        var_name = k.split("ref_")[1]
        assert allclose(l[var_name], m[k], v), (var_name, m[k] - l[var_name])

if __name__ == "__main__":
    test_regression()
