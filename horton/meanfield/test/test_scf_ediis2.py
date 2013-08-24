# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


import numpy as np
from horton import *
from horton.meanfield.test.common import check_scf_hf_cs_hf, check_scf_water_cs_hfs


def test_scf_ediis2_cs_hf():
    check_scf_hf_cs_hf(SCFWrapper('ediis2', threshold=1e-10, nvector=20))


def test_scf_ediis2_cs_hf_oda2():
    check_scf_hf_cs_hf(SCFWrapper('ediis2', threshold=1e-10, nvector=20, scf_step='oda2'))


def test_scf_ediis2_cs_hf_oda3():
    check_scf_hf_cs_hf(SCFWrapper('ediis2', threshold=1e-10, nvector=20, scf_step='oda3'))


def test_scf_ediis2_cs_hfs():
    check_scf_water_cs_hfs(SCFWrapper('ediis2', threshold=1e-10, nvector=20))


def test_scf_ediis2_cs_hfs_oda2():
    check_scf_water_cs_hfs(SCFWrapper('ediis2', threshold=1e-10, nvector=20, scf_step='oda3'))


def test_scf_ediis2_cs_hfs_oda3():
    check_scf_water_cs_hfs(SCFWrapper('ediis2', threshold=1e-10, nvector=20, scf_step='oda3'))
