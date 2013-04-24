# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
from horton.part.test.common import get_proatomdb_hf_sto3g, get_proatomdb_hf_lan


def check_names(names, wpart):
    for name in names:
        assert wpart.cache.has(name)


def check_water_hf_sto3g(scheme, local, expecting, **kwargs):
    proatomdb = get_proatomdb_hf_sto3g()

    # Get the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.update_dm('alpha')

    # Create a grid for the partitioning
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # Do the partitioning, both with local and global grids
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
    WPartClass = wpart_schemes[scheme]
    wpart = WPartClass(sys, grid, proatomdb, local, **kwargs)
    names = wpart.do_all()
    check_names(names, wpart)
    assert abs(wpart['charges'] - expecting).max() < 2e-3


def test_hirshfeld_water_hf_sto3g_local():
    expecting = np.array([-0.246171541212, 0.123092011074, 0.123079530138]) # from HiPart
    check_water_hf_sto3g('h', True, expecting)


def test_hirshfeld_water_hf_sto3g_global():
    expecting = np.array([-0.246171541212, 0.123092011074, 0.123079530138]) # from HiPart
    check_water_hf_sto3g('h', False, expecting)


def test_hirshfeld_i_water_hf_sto3g_local():
    expecting = np.array([-0.4214, 0.2107, 0.2107]) # From HiPart
    check_water_hf_sto3g('hi', True, expecting)


def test_hirshfeld_i_water_hf_sto3g_global():
    expecting = np.array([-0.4214, 0.2107, 0.2107]) # From HiPart
    check_water_hf_sto3g('hi', False, expecting)


def test_hirshfeld_e_water_hf_sto3g_local():
    expecting = np.array([-0.422794483125, 0.211390419810, 0.211404063315]) # From HiPart
    check_water_hf_sto3g('he', True, expecting)


def test_hirshfeld_e_water_hf_sto3g_global():
    expecting = np.array([-0.422794483125, 0.211390419810, 0.211404063315]) # From HiPart
    check_water_hf_sto3g('he', False, expecting)


def check_msa_hf_lan(scheme, local, expecting, **kwargs):
    proatomdb = get_proatomdb_hf_lan()

    # Get the molecule
    fn_fchk = context.get_fn('test/monosilicic_acid_hf_lan.fchk')
    sys = System.from_file(fn_fchk)

    # Create a grid for the partitioning
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # Do the partitioning, both with local and global grids
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
    WPartClass = wpart_schemes[scheme]
    wpart = WPartClass(sys, grid, proatomdb, local, **kwargs)
    wpart.do_charges()
    assert abs(wpart['charges'] - expecting).max() < 3e-3


def test_hirshfeld_msa_hf_lan_local():
    expecting = np.array([0.56175431, -0.30002709, -0.28602105, -0.28335086, -0.26832878,  0.13681904,  0.14535691,  0.14206876,  0.15097682])
    check_msa_hf_lan('h', True, expecting)


def test_hirshfeld_msa_hf_lan_global():
    expecting = np.array([0.56175431, -0.30002709, -0.28602105, -0.28335086, -0.26832878,  0.13681904,  0.14535691,  0.14206876,  0.15097682])
    check_msa_hf_lan('h', False, expecting)


def test_hirshfeld_i_msa_hf_lan_local():
    expecting = np.array([1.14305602, -0.52958298, -0.51787452, -0.51302759, -0.50033981, 0.21958586, 0.23189187, 0.22657354, 0.23938904])
    check_msa_hf_lan('hi', True, expecting)


def test_hirshfeld_i_msa_hf_lan_global():
    expecting = np.array([1.14305602, -0.52958298, -0.51787452, -0.51302759, -0.50033981, 0.21958586, 0.23189187, 0.22657354, 0.23938904])
    check_msa_hf_lan('hi', False, expecting)


def test_hirshfeld_e_msa_hf_lan_local():
    expecting = np.array([1.06135407, -0.51795437, -0.50626239, -0.50136175, -0.48867641, 0.22835963, 0.240736, 0.23528162, 0.24816043])
    check_msa_hf_lan('he', True, expecting)


def test_hirshfeld_e_msa_hf_lan_global():
    expecting = np.array([1.06135407, -0.51795437, -0.50626239, -0.50136175, -0.48867641, 0.22835963, 0.240736, 0.23528162, 0.24816043])
    check_msa_hf_lan('he', False, expecting)
