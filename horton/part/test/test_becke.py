#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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


from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_becke_n2_hfs_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 110), random_rotate=False, mode='only')
    dm_full = mol.get_dm_full()
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid.points)
    bp = BeckeWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens)
    bp.do_populations()
    assert abs(bp['populations'] - 7).max() < 1e-4
    bp.do_charges()
    assert abs(bp['charges']).max() < 1e-4
    bp.clear()
    with assert_raises(KeyError):
        bp['charges']
    bp.do_charges()
    assert abs(bp['populations'] - 7).max() < 1e-4
    assert abs(bp['charges']).max() < 1e-4


def test_becke_nonlocal_lih_hf_321g():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = IOData.from_file(fn_fchk)
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialGrid(rtf)
    dm_full = mol.get_dm_full()

    grid1 = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 110), random_rotate=False, mode='only')
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid1.points)
    bp1 = BeckeWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid1, moldens)

    grid2 = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, (rgrid, 110), random_rotate=False, mode='discard')
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid2.points)
    bp2 = BeckeWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid2, moldens, local=False)

    bp1.do_charges()
    bp2.do_charges()
    assert abs(bp1['charges'] - bp2['charges']).max() < 5e-4


def check_becke_azirine(key, expected):
    fn_fchk = context.get_fn('test/2h-azirine-%s.fchk' % key)
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False, mode='only')
    dm_full = mol.get_dm_full()
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid.points)
    bp = BeckeWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens)
    bp.do_charges()
    c = bp['charges']
    assert abs(c[0] - expected[0]) < 1e-3
    assert abs(c[2] - expected[1]) < 1e-3
    assert abs(c[5] - expected[2]) < 1e-3


def test_becke_azirine_ccd():
    check_becke_azirine('cc', [-0.0656538087277, -0.0770555290299, 0.123503410725])


def test_becke_azirine_cis():
    check_becke_azirine('ci', [-0.122893896731, -0.266685240737, 0.137147967309])


def test_becke_azirine_mp2():
    check_becke_azirine('mp2', [-0.0656579068849, -0.0761190062373, 0.126890127581])


def test_becke_azirine_mp3():
    check_becke_azirine('mp3', [-0.0665919182085, -0.0769654765789, 0.125587673579])


def test_becke_ch3_hf_sto3g():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False, mode='only')
    dm_full = mol.get_dm_full()
    dm_spin = mol.get_dm_spin()
    moldens = mol.obasis.compute_grid_density_dm(dm_full, grid.points)
    spindens = mol.obasis.compute_grid_density_dm(dm_spin, grid.points)
    bp = BeckeWPart(mol.coordinates, mol.numbers, mol.pseudo_numbers, grid, moldens, spindens)
    bp.do_all()
    sc = bp['spin_charges']
    assert abs(sc - [1.08458698, -0.02813376, -0.02813376, -0.02815979]).max() < 1e-3
