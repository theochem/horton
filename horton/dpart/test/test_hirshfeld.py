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


def get_proatomdb(qmin, qmax):
    # TODO: Make and use general-purpose proatomdb (in data directory)
    int1d = TrapezoidIntegrator1D()
    rtf = LogRTransform(1e-3, 1e1, 100)
    atgrid = AtomicGrid(np.zeros(3, float), rtf, int1d, 110, random_rotate=False, keep_subgrids=1)
    proatomdb = ProAtomDB.from_scratch([HartreeFock()], 'sto-3g', atgrid, [1,8], qmin=qmin, qmax=qmax)
    return proatomdb


def test_hirshfeld_water_hf_sto3g():
    proatomdb = get_proatomdb(0, 0)
    # Compute the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)

    # Create a grid for the partitionign
    int1d = TrapezoidIntegrator1D()
    rtf = LogRTransform(5e-4, 2e1, 120)
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=1)

    # do the partitioning
    hdp = HirshfeldDPart(grid, proatomdb)
    hdp.do_charges()
    expecting = np.array([-0.246, 0.123, 0.123]) # From HiPart
    assert abs(hdp['charges'] - expecting).max() < 1e-3


def test_hirshfeld_i_water_hf_sto3g():
    proatomdb = get_proatomdb(-1, 1)
    # Compute the molecule
    fn_fchk = context.get_fn('test/water_sto3g_hf_g03.fchk')
    sys = System.from_file(fn_fchk)

    # Create a grid for the partitionign
    int1d = TrapezoidIntegrator1D()
    rtf = LogRTransform(5e-4, 2e1, 120)
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False, keep_subgrids=1)

    # do the partitioning
    hdp = HirshfeldIDPart(grid, proatomdb, 1e-4)
    hdp.do_charges()
    expecting = np.array([-0.4214, 0.2107, 0.2107]) # From HiPart
    assert abs(hdp['charges'] - expecting).max() < 1e-3
