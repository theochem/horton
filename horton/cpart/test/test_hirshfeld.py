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
from horton.cpart.test.common import get_fake_co, get_fake_pseudo_oo


def check_names(names, cpart):
    for name in names:
        assert cpart._cache.has(name)


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
    proatomdb = ProAtomDB.from_refatoms(atgrid, numbers=[8,14], max_kation=0, max_anion=0)

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_jbw_coarse.h5') as store:
        cpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb, store)
        names = cpart.do_all()
        check_names(names, cpart)
        wcor = cpart['wcor']
        assert abs(cpart['populations'].sum() - ui_grid.integrate(wcor, mol_dens)) < 1e-10


def test_hirshfeld_fake():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_fake.h5') as store:
        cpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb, store)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(abs(charges).mean() - 0.114) < 1e-3


def test_hirshfeld_fake_pseudo():
    sys, ui_grid, mol_dens, proatomdb = get_fake_pseudo_oo()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_fake_pseudo.h5') as store:
        cpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb, store)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(charges[0] - 0.2119886) < 1e-3


def test_hirshfeld_fake_pseudo_smooth():
    sys, ui_grid, mol_dens, proatomdb = get_fake_pseudo_oo()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_fake_pseudo_smooth.h5') as store:
        cpart = HirshfeldCPart(sys, ui_grid, mol_dens, proatomdb, store, smooth=True)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(charges[0] - 0.2119886) < 1e-3


def check_proatom_splines(cpart):
    for index in xrange(cpart.system.natom):
        center = cpart.system.coordinates[index]
        spline = cpart.get_proatom_spline(index)
        array1 = cpart.ui_grid.zeros()
        cpart.ui_grid.eval_spline(spline, center, array1)
        array2 = cpart.ui_grid.zeros()
        cpart.compute_proatom(index, array2)
        assert abs(array1).max() != 0.0
        assert abs(array1 - array2).max() < 1e-5


def test_hirshfeld_i_fake():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_i_fake') as store:
        cpart = HirshfeldICPart(sys, ui_grid, mol_dens, proatomdb, store)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(abs(charges).mean() - 0.438) < 1e-3
        check_proatom_splines(cpart)


def test_hirshfeld_i_fake_pseudo():
    sys, ui_grid, mol_dens, proatomdb = get_fake_pseudo_oo()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_i_fake_pseudo') as store:
        cpart = HirshfeldICPart(sys, ui_grid, mol_dens, proatomdb, store)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(charges[0] - 0.40262645) < 1e-3
        check_proatom_splines(cpart)


def test_hirshfeld_i_fake_pseudo_smooth():
    sys, ui_grid, mol_dens, proatomdb = get_fake_pseudo_oo()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_i_fake_pseudo_smooth') as store:
        cpart = HirshfeldICPart(sys, ui_grid, mol_dens, proatomdb, store, smooth=True)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(charges[0] - 0.40262645) < 1e-3
        check_proatom_splines(cpart)


def test_hirshfeld_e_fake():
    sys, ui_grid, mol_dens, proatomdb = get_fake_co()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_e_fake') as store:
        cpart = HirshfeldECPart(sys, ui_grid, mol_dens, proatomdb, store)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(abs(charges).mean() - 0.419) < 1e-3
        check_proatom_splines(cpart)


def test_hirshfeld_e_fake_pseudo():
    sys, ui_grid, mol_dens, proatomdb = get_fake_pseudo_oo()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_e_fake_pseudo') as store:
        cpart = HirshfeldECPart(sys, ui_grid, mol_dens, proatomdb, store)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(charges[0] - 0.400) < 1e-3
        check_proatom_splines(cpart)


def test_hirshfeld_e_fake_pseudo_smooth():
    sys, ui_grid, mol_dens, proatomdb = get_fake_pseudo_oo()

    # Run the partitioning
    with ArrayStore.from_mode('core', 'horton.test.test_hirshfeld.test_hirshfeld_e_fake_pseudo_smooth') as store:
        cpart = HirshfeldECPart(sys, ui_grid, mol_dens, proatomdb, store, smooth=True)
        cpart.do_charges()
        charges = cpart['charges']
        assert abs(charges.sum()) < 1e-2
        assert abs(charges[0] - 0.400) < 1e-3
        check_proatom_splines(cpart)
