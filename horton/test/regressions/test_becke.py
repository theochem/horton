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

from numpy import array, allclose
from nose.plugins.attrib import attr

from horton import context


@attr('regression_check')
def test_becke():
    ref_result_charges = array([-0.13537344,  0.06766372,  0.06767448])

    results = ['ref_result_charges']
    thresholds = {'ref_result_charges': 1e-08}

    test_path = context.get_fn("examples/wpart/becke.py")
    with open(test_path) as fh:
        exec fh

    l = locals()
    for r in results:
        var_name = r.split("ref_")[-1]
        assert allclose(l[var_name], l[r], thresholds[r]), l[r] - l[var_name]