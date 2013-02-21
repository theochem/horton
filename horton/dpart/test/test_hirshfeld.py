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
from horton.dpart.test.common import get_proatomdb_hf_sto3g, get_proatomdb_hf_lan


# TODO: reduce duplicate code

def test_hirshfeld_water_hf_sto3g():
    proatomdb = get_proatomdb_hf_sto3g()

    # Get the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.update_dm('alpha')

    # Create a grid for the partitioning
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # Do the partitioning, both with local and global grids
    for local in True, False:
        grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
        dpart = HirshfeldDPart(grid, proatomdb, local)
        dpart.do_charges()
        expecting = np.array([-0.246171541212, 0.123092011074, 0.123079530138]) # from HiPart
        assert abs(dpart['charges'] - expecting).max() < 2e-3


def test_hirshfeld_msa_hf_lan():
    proatomdb = get_proatomdb_hf_lan()

    # Get the molecule
    fn_fchk = context.get_fn('test/monosilicic_acid_hf_lan.fchk')
    sys = System.from_file(fn_fchk)

    # Create a grid for the partitioning
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # Do the partitioning, both with local and global grids
    for local in True, False:
        grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
        dpart = HirshfeldDPart(grid, proatomdb, local)
        dpart.do_charges()
        expecting = np.array([0.56175431, -0.30002709, -0.28602105, -0.28335086, -0.26832878,  0.13681904,  0.14535691,  0.14206876,  0.15097682])
        assert abs(dpart['charges'] - expecting).max() < 2e-3


def test_hirshfeld_i_water_hf_sto3g():
    proatomdb = get_proatomdb_hf_sto3g()

    # Get the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.update_dm('alpha')

    # Create a grid for the partitioning
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # Do the partitioning, both with local and global grids
    for local in True, False:
        grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
        dpart = HirshfeldIDPart(grid, proatomdb, local, 1e-4)
        dpart.do_charges()
        expecting = np.array([-0.4214, 0.2107, 0.2107]) # From HiPart
        assert abs(dpart['charges'] - expecting).max() < 1e-3


def test_hirshfeld_i_msa_hf_lan():
    proatomdb = get_proatomdb_hf_lan()

    # Get the molecule
    fn_fchk = context.get_fn('test/monosilicic_acid_hf_lan.fchk')
    sys = System.from_file(fn_fchk)

    # Create a grid for the partitioning
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # Do the partitioning, both with local and global grids
    for local in True, False:
        grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
        dpart = HirshfeldIDPart(grid, proatomdb, local)
        dpart.do_charges()
        expecting = np.array([1.14305602, -0.52958298, -0.51787452, -0.51302759, -0.50033981, 0.21958586, 0.23189187, 0.22657354, 0.23938904])
        assert abs(dpart['charges'] - expecting).max() < 2e-3


def test_hirshfeld_e_water_hf_sto3g():
    proatomdb = get_proatomdb_hf_sto3g()

    # Get the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.update_dm('alpha')

    # Create a grid for the partitioning
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # Do the partitioning, both with local and global grids
    for local in True, False:
        grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
        dpart = HirshfeldEDPart(grid, proatomdb, local, 1e-4)
        dpart.do_charges()
        expecting = np.array([-0.422794483125, 0.211390419810, 0.211404063315]) # From HiPart
        assert abs(dpart['charges'] - expecting).max() < 1e-3
