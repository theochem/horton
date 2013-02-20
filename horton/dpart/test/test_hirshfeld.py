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
from horton.dpart.test.common import get_proatomdb_ref


def test_hirshfeld_water_hf_sto3g():
    proatomdb = get_proatomdb_ref([1, 8], 0, 0)
    # Compute the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.update_dm('alpha')

    # Create a grid for the partitionign
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # do the partitioning, both with local and global grids
    for local in True, False:
        grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
        hdp = HirshfeldDPart(grid, proatomdb, local)
        hdp.do_charges()
        expecting = np.array([-0.509,  0.253,  0.254]) # From Horton
        assert abs(hdp['charges'] - expecting).max() < 2e-3


def test_hirshfeld_i_water_hf_sto3g():
    # TODO: find other test case or other pro-atoms
    from nose.plugins.skip import SkipTest
    raise SkipTest

    proatomdb = get_proatomdb_ref([1, 8], 2, 2)
    # Compute the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.update_dm('alpha')

    # Create a grid for the partitionign
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # do the partitioning
    for local in True, False:
        grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
        hdp = HirshfeldIDPart(grid, proatomdb, local, 1e-4)
        hdp.do_charges()
        expecting = np.array([-0.4214, 0.2107, 0.2107]) # From HiPart
        assert abs(hdp['charges'] - expecting).max() < 1e-3


def test_hirshfeld_e_water_hf_sto3g():
    # TODO: find other test case or other pro-atoms
    from nose.plugins.skip import SkipTest
    raise SkipTest

    proatomdb = get_proatomdb_ref([1, 8], 1, 1)
    # Compute the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.update_dm('alpha')

    # Create a grid for the partitionign
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(5e-4, 2e1, 120)

    # do the partitioning
    for local in True, False:
        grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=int(local))
        hdp = HirshfeldEDPart(grid, proatomdb, local, 1e-4)
        hdp.do_charges()
        expecting = np.array([-0.422794483125, 0.211390419810, 0.211404063315]) # From HiPart
        assert abs(hdp['charges'] - expecting).max() < 1e-3
