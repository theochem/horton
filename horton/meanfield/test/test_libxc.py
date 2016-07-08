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


import numpy as np
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_interpolation


def test_cubic_interpolation_c_pbe_cs():
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCGGA('c_pbe'),
        ]),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha])


def test_cubic_interpolation_x_pbe_cs():
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCGGA('x_pbe'),
        ]),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha])


def test_cubic_interpolation_hfs_cs():
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCLDA('x'),
        ]),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha])


def test_cubic_interpolation_x_pbe_c_vwn_cs():
    # mixing of GGA and LDA
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCGGA('x_pbe'),
            RLibXCLDA('c_vwn'),
        ]),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha])


def test_cubic_interpolation_c_vwn_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        RGridGroup(mol.obasis, grid, [
            RLibXCLDA('c_vwn'),
        ]),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha])


def test_cubic_interpolation_o3lyp_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    libxc_term = RLibXCHybridGGA('xc_o3lyp')
    terms = [
        RGridGroup(mol.obasis, grid, [libxc_term]),
        RExchangeTerm(er, 'x_hf', libxc_term.get_exx_fraction()),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha])


def test_cubic_interpolation_x_tpss_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        RGridGroup(mol.obasis, grid, [RLibXCMGGA('x_tpss')]),
    ]
    ham = REffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha])


def test_cubic_interpolation_c_pbe_os():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        UGridGroup(mol.obasis, grid, [
            ULibXCGGA('c_pbe'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha, mol.exp_beta])


def test_cubic_interpolation_x_pbe_os():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        UGridGroup(mol.obasis, grid, [
            ULibXCGGA('x_pbe'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha, mol.exp_beta])


def test_cubic_interpolation_hfs_os():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        UGridGroup(mol.obasis, grid, [
            ULibXCLDA('x'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha, mol.exp_beta])


def test_cubic_interpolation_x_pbe_c_vwn_os():
    # mixing of LDA and GGA
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        UGridGroup(mol.obasis, grid, [
            ULibXCGGA('x_pbe'),
            ULibXCLDA('c_vwn'),
        ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha, mol.exp_beta])


def test_cubic_interpolation_o3lyp_os():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    libxc_term = ULibXCHybridGGA('xc_o3lyp')
    terms = [
        UGridGroup(mol.obasis, grid, [libxc_term]),
        UExchangeTerm(er, 'x_hf', libxc_term.get_exx_fraction()),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha, mol.exp_beta])


def test_cubic_interpolation_x_tpss_os():
    # mixing of LDA and GGA
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    terms = [
        UGridGroup(mol.obasis, grid, [ULibXCMGGA('x_tpss'), ]),
    ]
    ham = UEffHam(terms)
    check_interpolation(ham, mol.lf, olp, kin, na, [mol.exp_alpha, mol.exp_beta])


def test_hyb_gga_exx_fraction():
    # xc_pbeh = The PBE0 functional
    t1 = RLibXCHybridGGA('xc_pbeh')
    assert t1.get_exx_fraction() == 0.25
    t2 = ULibXCHybridGGA('xc_pbeh')
    assert t2.get_exx_fraction() == 0.25


def test_lda_c_vwn_present():
    t1 = RLibXCLDA('c_vwn')     # The VWN 5 functional
    assert t1._libxc_wrapper.key == 'lda_c_vwn'
    t2 = RLibXCLDA('c_vwn_4')   # The VWN 4 functional
    assert t2._libxc_wrapper.key == 'lda_c_vwn_4'


def test_info():
    t = RLibXCWrapper('lda_x')
    assert t.key == 'lda_x'
    assert t.name is not None
    assert t.number is not None
    assert t.kind is not None
    assert t.family is not None
    assert t.refs is not None


def test_hyb_mgga_exx_fraction():
    # xc_tpssh = The TPSS functional with exact exchange
    t1 = RLibXCHybridMGGA('xc_tpssh')
    assert t1.get_exx_fraction() == 0.1
    t2 = ULibXCHybridMGGA('xc_tpssh')
    assert t2.get_exx_fraction() == 0.1
