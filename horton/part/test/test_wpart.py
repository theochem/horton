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


import numpy as np
from nose.plugins.attrib import attr

from horton.grid import ExpRTransform, RadialGrid, BeckeMolGrid
from .. utils import wpart_schemes
from .common import (load_molecule_npz, check_names, check_proatom_splines,
                     get_proatomdb_hf_sto3g, get_proatomdb_hf_lan)


def check_water_hf_sto3g(scheme, expecting, needs_padb=True, **kwargs):
    if needs_padb:
        proatomdb = get_proatomdb_hf_sto3g()
        kwargs['proatomdb'] = proatomdb

    # load molecule data
    coords, nums, pseudo_nums, dens, points = load_molecule_npz('water_sto3g_hf_g03_fchk_exp:5e-4:2e1:120:110.npz')
    # Create a grid for the partitioning
    rtf = ExpRTransform(5e-4, 2e1, 120)
    rgrid = RadialGrid(rtf)
    mode = 'only' if kwargs.get('local', True) else 'discard'
    grid = BeckeMolGrid(coords, nums, pseudo_nums, (rgrid, 110), random_rotate=False, mode=mode)
    # check the grid points against stored points on which density is evaluated
    assert (abs(points - grid.points) < 1.e-6).all()
    # Do the partitioning
    WPartClass = wpart_schemes(scheme)
    wpart = WPartClass(coords, nums, pseudo_nums, grid, dens,  **kwargs)
    names = wpart.do_all()
    check_names(names, wpart)
    assert abs(wpart['charges'] - expecting).max() < 2e-3
    assert abs(wpart['charges'] - wpart['cartesian_multipoles'][:,0]).max() < 1e-3
    assert abs(wpart['charges'] - wpart['pure_multipoles'][:,0]).max() < 1e-3

    check_proatom_splines(wpart)
    return wpart


def test_hirshfeld_water_hf_sto3g_local():
    expecting = np.array([-0.246171541212, 0.123092011074, 0.123079530138]) # from HiPart
    check_water_hf_sto3g('h', expecting, local=True)


def test_hirshfeld_water_hf_sto3g_global():
    expecting = np.array([-0.246171541212, 0.123092011074, 0.123079530138]) # from HiPart
    check_water_hf_sto3g('h', expecting, local=False)


def test_hirshfeld_i_water_hf_sto3g_local():
    expecting = np.array([-0.4214, 0.2107, 0.2107]) # From HiPart
    check_water_hf_sto3g('hi', expecting, local=True)


def test_hirshfeld_i_water_hf_sto3g_global():
    expecting = np.array([-0.4214, 0.2107, 0.2107]) # From HiPart
    check_water_hf_sto3g('hi', expecting, local=False)


def test_is_water_hf_sto3g():
    expecting = np.array([-0.490017586929, 0.245018706885, 0.244998880045]) # From HiPart
    check_water_hf_sto3g('is', expecting, needs_padb=False)


def test_mbis_water_hf_sto3g():
    expecting = np.array([-0.61891067, 0.3095756, 0.30932584])
    wpart = check_water_hf_sto3g('mbis', expecting, needs_padb=False)
    assert (wpart['charges'] == wpart['valence_charges'] + wpart['core_charges']).all()
    assert (wpart['core_charges'] > 0).all()
    assert (wpart['valence_charges'] < 0).all()
    assert (wpart['valence_widths'] > 0).all()


def check_msa_hf_lan(scheme, expecting, needs_padb=True, **kwargs):
    if needs_padb:
        proatomdb = get_proatomdb_hf_lan()
        kwargs['proatomdb'] = proatomdb

    # load molecule data
    coords, nums, pseudo_nums, dens, points = load_molecule_npz('monosilicic_acid_hf_lan_fchk_exp:5e-4:2e1:120:110.npz')
    # Create a grid for the partitioning
    rtf = ExpRTransform(5e-4, 2e1, 120)
    rgrid = RadialGrid(rtf)
    # Do the partitioning, both with local and global grids
    mode = 'only' if kwargs.get('local', True) else 'discard'
    grid = BeckeMolGrid(coords, nums, pseudo_nums, (rgrid, 110), random_rotate=False, mode=mode)
    # check the grid points against stored points on which density is evaluated
    assert (abs(points - grid.points) < 1.e-6).all()
    # Do the partitioning
    WPartClass = wpart_schemes(scheme)
    wpart = WPartClass(coords, nums, pseudo_nums, grid, dens, **kwargs)
    wpart.do_charges()
    assert abs(wpart['charges'] - expecting).max() < 4e-3

    check_proatom_splines(wpart)


def test_hirshfeld_msa_hf_lan_local():
    expecting = np.array([0.56175431, -0.30002709, -0.28602105, -0.28335086, -0.26832878,  0.13681904,  0.14535691,  0.14206876,  0.15097682])
    check_msa_hf_lan('h', expecting, local=True)


@attr('slow')
def test_hirshfeld_msa_hf_lan_global():
    expecting = np.array([0.56175431, -0.30002709, -0.28602105, -0.28335086, -0.26832878,  0.13681904,  0.14535691,  0.14206876,  0.15097682])
    check_msa_hf_lan('h', expecting, local=False)


@attr('slow')
def test_hirshfeld_i_msa_hf_lan_local():
    expecting = np.array([1.14305602, -0.52958298, -0.51787452, -0.51302759, -0.50033981, 0.21958586, 0.23189187, 0.22657354, 0.23938904])
    check_msa_hf_lan('hi', expecting, local=True)


@attr('slow')
def test_hirshfeld_i_msa_hf_lan_global():
    expecting = np.array([1.14305602, -0.52958298, -0.51787452, -0.51302759, -0.50033981, 0.21958586, 0.23189187, 0.22657354, 0.23938904])
    check_msa_hf_lan('hi', expecting, local=False)


@attr('slow')
def test_is_msa_hf_lan():
    expecting = np.array([1.1721364, -0.5799622, -0.5654549, -0.5599638, -0.5444145, 0.2606699, 0.2721848, 0.2664377, 0.2783666]) # from HiPart
    check_msa_hf_lan('is', expecting, needs_padb=False)
