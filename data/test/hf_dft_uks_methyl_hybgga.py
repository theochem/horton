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
    ref_result_kin = 39.34908273259613
    ref_result_hartree = 28.101087291309366
    ref_result_energy = -39.829288592008702
    ref_result_ex = -1.2302308724933995
    ref_result_grid = -5.206316404353287
    ref_result_exp_alpha = array([-10.20251592,  -0.66406422,  -0.40274032,  -0.4027181 ,
        -0.22486684,   0.0888693 ,   0.16219268,   0.16220226,
         0.5332166 ,   0.56768584,   0.56769077,   0.68331588,
         0.86344917,   0.86351231,   0.92072183,   1.65654082,
         1.65655893,   1.95383736,   2.10939831,   2.10951572])
    ref_result_exp_beta = array([-10.18784099,  -0.62229874,  -0.39383725,  -0.39381147,
        -0.04713796,   0.10635446,   0.16868605,   0.16869399,
         0.57658235,   0.57659161,   0.6121897 ,   0.70331622,
         0.86667792,   0.86674732,   0.95060174,   1.72983425,
         1.72985687,   2.04447856,   2.12113476,   2.12123425])
    ref_result_nn = 9.0797849426636361
    ref_result_ne = -109.92269628173115

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
