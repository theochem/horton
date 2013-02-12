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


def test_hirshfeld_jbw_coarse():
    # This test is not supposed to generate meaningful numbers. The cube data
    # is too coarse and the reference atoms may have little similarities with
    # the DFT density.

    # Load the cube file
    fn_cube = context.get_fn('test/jbw_coarse_aedens.cube')
    sys = System.from_file(fn_cube)
    mol_dens = sys.props['cube_data']
    ui_grid = sys.props['ui_grid']

    # Load some pro-atoms
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rtf, int1d, 110), keep_subgrids=1)
    proatomdb = ProAtomDB.from_refatoms(atgrid, numbers=[8,14], qmax=0)

    # Run the partitioning
    hcpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb)
    hcpart.do_charges()
    assert abs(hcpart['charges'].sum() + ui_grid.integrate(hcpart['rel_mol_dens'])) < 1e-10


def get_fake_co():
    # Define system
    numbers = np.array([6, 8])
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.132]])
    sys = System(coordinates, numbers)

    # Load some pro-atoms
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rtf, int1d, 110), keep_subgrids=1)
    proatomdb = ProAtomDB.from_refatoms(atgrid, numbers=[6,8], qmax=1)

    # Make fake cube data
    origin = np.array([-3.0, -3.0, -3.0])
    rvecs = np.identity(3, float)*0.1
    shape = np.array([60, 60, 60+22])
    ui_grid = UniformIntGrid(origin, Cell(rvecs), shape)

    mol_dens = np.zeros(ui_grid.shape)
    tmp = np.zeros(ui_grid.shape)
    setup = [
        (0, 5, 0.5),
        (0, 6, 0.4),
        (0, 7, 0.1),
        (1, 7, 0.1),
        (1, 8, 0.4),
        (1, 9, 0.5),
    ]
    for i, pop, frac in setup:
        n = sys.numbers[i]
        c = sys.coordinates[i]
        # TODO: can be made more efficient by scaling the spline function
        spline = proatomdb.get_hirshfeld_i_proatom_fn(n, pop)
        tmp[:] = 0.0
        ui_grid.eval_spline(spline, c, tmp)
        tmp *= frac
        mol_dens += tmp

    return sys, ui_grid, mol_dens, proatomdb


def test_hirshfeld_fake():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    # Run the partitioning
    hicpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb)
    hicpart.do_charges()
    charges = hicpart['charges']
    assert charges.sum() < 1e-3
    assert abs(charges[0] - 0.11229602) < 1e-3


def test_hirshfeld_i_fake():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    # Run the partitioning
    hicpart = HirshfeldICPart(sys, ui_grid, mol_dens, proatomdb)
    hicpart.do_charges()
    charges = hicpart['charges']
    assert charges.sum() < 1e-3
    assert abs(charges[0] - 0.433) < 1e-3
