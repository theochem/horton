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
from horton.cpart.test.common import get_fake_co


def test_cpart_hirshfeld_jbw_coarse():
    # This test is not supposed to generate meaningful numbers. The cube data
    # is too coarse and the reference atoms may have little similarities with
    # the DFT density.

    # Load the cube file
    fn_cube = context.get_fn('test/jbw_coarse_aedens.cube')
    sys = System.from_file(fn_cube)
    del sys.props['pseudo_numbers']
    mol_dens = sys.props['cube_data']
    ui_grid = sys.props['ui_grid']

    # Load some pro-atoms
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rtf, int1d, 110), keep_subgrids=1)
    proatomdb = ProAtomDB.from_refatoms(atgrid, numbers=[8,14], qmax=0)

    # Run the partitioning
    cpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb, False)
    cpart.do_charges()
    cpart.do_dipoles()
    cpart.do_volumes()
    assert abs(cpart['populations'].sum() - cpart._integrate(cpart['moldens'])) < 1e-10
    assert cpart['volumes'].min() > 0


def test_cpart_hirshfeld_fake():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    # Run the partitioning
    cpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb, False)
    cpart.do_charges()
    charges = cpart['charges']
    assert charges.sum() < 1e-3
    assert abs(charges[0] - 0.112) < 1e-3


def test_cpart_hirshfeld_fake_smooth():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co(smooth=True)

    # Run the partitioning
    cpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb, True)
    cpart.do_charges()
    charges = cpart['charges']
    assert charges.sum() < 1e-3
    assert abs(charges[0] - 0.120) < 1e-3


def test_cpart_hirshfeld_i_fake():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    # Run the partitioning
    cpart = HirshfeldICPart(sys, ui_grid, mol_dens, proatomdb, False)
    cpart.do_charges()
    charges = cpart['charges']
    assert charges.sum() < 1e-3
    assert abs(charges[0] - 0.431) < 1e-3


def test_cpart_hirshfeld_i_fake_smooth():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co(smooth=True)

    # Run the partitioning
    cpart = HirshfeldICPart(sys, ui_grid, mol_dens, proatomdb, True)
    cpart.do_charges()
    charges = cpart['charges']
    assert charges.sum() < 1e-3
    assert abs(charges[0] - 0.113) < 1e-3


def test_cpart_hirshfeld_e_fake():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    # Run the partitioning
    cpart = HirshfeldECPart(sys, ui_grid, mol_dens, proatomdb, False)
    cpart.do_charges()
    charges = cpart['charges']
    assert charges.sum() < 1e-3
    assert abs(charges[0] - 0.393) < 1e-3
