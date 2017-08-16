#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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

from .common import load_molecule_npz, get_fn
from .. becke import BeckeWPart
from horton.grid import ExpRTransform, RadialGrid, BeckeMolGrid


def test_becke_n2_hfs_sto3g():
    # load molecule data
    coords, nums, pseudo_nums, dens, points = load_molecule_npz('n2_hfs_sto3g_fchk_exp:1e-3:1e1:100:110.npz')
    # make grid
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(coords, nums, pseudo_nums, (rgrid, 110), random_rotate=False, mode='only')
    # check the grid points against stored points on which density is evaluated
    assert (abs(points - grid.points) < 1.e-6).all()
    # Becke partitioning
    bp = BeckeWPart(coords, nums, pseudo_nums, grid, dens)
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
    # load molecule data
    coords, nums, pnums, dens, points = load_molecule_npz('li_h_3-21G_hf_g09_fchk_exp:1e-3:1e1:100:110.npz')

    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialGrid(rtf)

    grid1 = BeckeMolGrid(coords, nums, pnums, (rgrid, 110), random_rotate=False, mode='only')
    bp1 = BeckeWPart(coords, nums, pnums, grid1, dens)
    # check the grid points against stored points on which density is evaluated
    assert (abs(points - grid1.points) < 1.e-6).all()

    grid2 = BeckeMolGrid(coords, nums, pnums, (rgrid, 110), random_rotate=False, mode='discard')
    bp2 = BeckeWPart(coords, nums, pnums, grid2, dens, local=False)
    # check the grid points against stored points on which density is evaluated
    assert (abs(points - grid2.points) < 1.e-6).all()

    bp1.do_charges()
    bp2.do_charges()
    assert abs(bp1['charges'] - bp2['charges']).max() < 5e-4


def check_becke_azirine(key, expected):
    # load molecule data
    coords, nums, pseudo_nums, dens, points = load_molecule_npz('2h-azirine-%s-fchk-medium.npz' % key)
    # make grid
    grid = BeckeMolGrid(coords, nums, pseudo_nums, random_rotate=False, mode='only', agspec='medium')
    # check the grid points against stored points on which density is evaluated
    assert (abs(points - grid.points) < 1.e-6).all()
    # Becke partitioning
    bp = BeckeWPart(coords, nums, pseudo_nums, grid, dens)
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
    # load molecule data
    coords, nums, pnums, dens, spindens, points = load_molecule_npz('ch3_hf_sto3g_fchk_medium.npz', True)
    grid = BeckeMolGrid(coords, nums, pnums, random_rotate=False, mode='only', agspec='medium')
    bp = BeckeWPart(coords, nums, pnums, grid, dens, spindens)
    bp.do_all()
    sc = bp['spin_charges']
    assert abs(sc - [1.08458698, -0.02813376, -0.02813376, -0.02815979]).max() < 1e-3
