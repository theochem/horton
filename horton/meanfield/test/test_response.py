# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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



import numpy as np

from horton.meanfield.test.common import get_obasis, load_orbs_alpha, load_orbs_beta
from horton.part import get_mulliken_operators
from .. import compute_noninteracting_response


def check_response(fname):
    operators = get_mulliken_operators(get_obasis(fname))
    orbs = [load_orbs_alpha(fname)]
    try:
        orbs.append(load_orbs_beta(fname))
    except IOError:
        pass

    for orb in orbs:
        response = compute_noninteracting_response(orb, operators)
        assert np.isfinite(response).all()
        assert abs(response.sum(axis=0)).max() < 1e-3
        evals = np.linalg.eigvalsh(response)
        assert (evals[:2] < 0).all()
        assert abs(evals[-1]) < 1e-10


def test_response_water_sto3g():
    fname = 'water_sto3g_hf_g03_fchk'
    check_response(fname)


def test_response_h3_321g():
    fname = 'h3_hfs_321g_fchk'
    check_response(fname)

# TODO: Move to higher level test
# def test_response_benzene():
#     fname = 'benzene_sto3g_fchk'
#     check_response(fname)
