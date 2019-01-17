# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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

from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_hf_cs_hf, check_lih_os_hf, \
    check_water_cs_hfs, check_n2_cs_hfs, check_h3_os_hfs, check_h3_os_pbe, \
    check_co_cs_pbe, check_vanadium_sc_hf, check_water_cs_m05, \
    check_methyl_os_tpss


def test_hf_cs_hf():
    check_hf_cs_hf(CDIISSCFSolver(threshold=1e-7))


def test_lih_os_hf():
    check_lih_os_hf(CDIISSCFSolver(threshold=1e-7))


def test_water_cs_hfs():
    check_water_cs_hfs(CDIISSCFSolver(threshold=1e-6))


def test_n2_cs_hfs():
    check_n2_cs_hfs(CDIISSCFSolver(threshold=1e-6))


def test_h3_os_hfs():
    check_h3_os_hfs(CDIISSCFSolver(threshold=1e-6))


def test_co_cs_pbe():
    check_co_cs_pbe(CDIISSCFSolver(threshold=1e-5))


def test_h3_os_pbe():
    check_h3_os_pbe(CDIISSCFSolver(threshold=1e-6))


def test_vanadium_sc_hf():
    with assert_raises(NoSCFConvergence):
        check_vanadium_sc_hf(CDIISSCFSolver(threshold=1e-10, maxiter=10))


def test_water_cs_m05():
    check_water_cs_m05(CDIISSCFSolver(threshold=1e-6))


def test_methyl_os_tpss():
    check_methyl_os_tpss(CDIISSCFSolver(threshold=1e-5))
