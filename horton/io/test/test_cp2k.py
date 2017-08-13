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
"""Tests for CP2K ATOM output reader."""

from nose.tools import assert_raises

from . common import get_fn, truncated_file, check_orthonormal

from .. iodata import IOData
from .. overlap import compute_overlap



# TODO: add more obasis tests?


def check_orthonormality(mol):
    """Helper function to test if the orbitals are orthonormal."""
    olp = compute_overlap(**mol.obasis)
    check_orthonormal(mol.orb_alpha_occs, mol.orb_alpha_coeffs, olp)
    if hasattr(mol, 'orb_beta'):
        check_orthonormal(mol.orb_beta_occs, mol.orb_beta_coeffs, olp)


def test_atom_si_uks():
    fn_out = get_fn('atom_si.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [14]).all()
    assert (mol.pseudo_numbers == [4]).all()
    assert (mol.orb_alpha_occs == [1, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]).all()
    assert (mol.orb_beta_occs == [1, 0, 0, 0]).all()
    assert abs(mol.orb_alpha_energies - [-0.398761, -0.154896, -0.154896, -0.154896]).max() < 1e-4
    assert abs(mol.orb_beta_energies - [-0.334567, -0.092237, -0.092237, -0.092237]).max() < 1e-4
    assert abs(mol.energy - -3.761587698067) < 1e-10
    assert (mol.obasis["shell_types"] == [0, 0, 1, 1, -2]).all()
    check_orthonormality(mol)


def test_atom_o_rks():
    fn_out = get_fn('atom_om2.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [8]).all()
    assert (mol.pseudo_numbers == [6]).all()
    assert (mol.orb_alpha_occs == [1, 1, 1, 1]).all()
    assert abs(mol.orb_alpha_energies - [0.102709, 0.606458, 0.606458, 0.606458]).max() < 1e-4
    assert not hasattr(mol, 'orb_beta')
    assert abs(mol.energy - -15.464982778766) < 1e-10
    assert (mol.obasis["shell_types"] == [0, 0, 1, 1, -2]).all()
    check_orthonormality(mol)


def test_carbon_gs_ae_contracted():
    fn_out = get_fn('carbon_gs_ae_contracted.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [6]).all()
    assert (mol.pseudo_numbers == [6]).all()
    assert (mol.orb_alpha_occs == [1, 1, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]).all()
    assert (mol.orb_alpha_energies == [-10.058194, -0.526244, -0.214978,
                                       -0.214978, -0.214978]).all()
    assert (mol.orb_beta_occs == [1, 1, 0, 0, 0]).all()
    assert (mol.orb_beta_energies == [-10.029898, -0.434300, -0.133323,
                                      -0.133323, -0.133323]).all()
    assert mol.energy == -37.836423363057
    check_orthonormality(mol)


def test_carbon_gs_ae_uncontracted():
    fn_out = get_fn('carbon_gs_ae_uncontracted.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [6]).all()
    assert (mol.pseudo_numbers == [6]).all()
    assert (mol.orb_alpha_occs == [1, 1, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]).all()
    assert (mol.orb_alpha_energies == [-10.050076, -0.528162, -0.217626,
                                       -0.217626, -0.217626]).all()
    assert (mol.orb_beta_occs == [1, 1, 0, 0, 0]).all()
    assert (mol.orb_beta_energies == [-10.022715, -0.436340, -0.137135,
                                      -0.137135, -0.137135]).all()
    assert mol.energy == -37.842552743398
    check_orthonormality(mol)


def test_carbon_gs_pp_contracted():
    fn_out = get_fn('carbon_gs_pp_contracted.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [6]).all()
    assert (mol.pseudo_numbers == [4]).all()
    assert (mol.orb_alpha_occs == [1, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]).all()
    assert (mol.orb_alpha_energies == [-0.528007, -0.219974, -0.219974, -0.219974]).all()
    assert (mol.orb_beta_occs == [1, 0, 0, 0]).all()
    assert (mol.orb_beta_energies == [-0.429657, -0.127060, -0.127060, -0.127060]).all()
    assert mol.energy == -5.399938535844
    check_orthonormality(mol)


def test_carbon_gs_pp_uncontracted():
    fn_out = get_fn('carbon_gs_pp_uncontracted.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [6]).all()
    assert (mol.pseudo_numbers == [4]).all()
    assert (mol.orb_alpha_occs == [1, 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]).all()
    assert (mol.orb_alpha_energies == [-0.528146, -0.219803, -0.219803, -0.219803]).all()
    assert (mol.orb_beta_occs == [1, 0, 0, 0]).all()
    assert (mol.orb_beta_energies == [-0.429358, -0.126411, -0.126411, -0.126411]).all()
    assert mol.energy == -5.402288849332
    check_orthonormality(mol)


def test_carbon_sc_ae_contracted():
    fn_out = get_fn('carbon_sc_ae_contracted.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [6]).all()
    assert (mol.pseudo_numbers == [6]).all()
    assert (mol.orb_alpha_occs == [1, 1, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]).all()
    assert (mol.orb_alpha_energies == [-10.067251, -0.495823, -0.187878,
                                       -0.187878, -0.187878]).all()
    assert not hasattr(mol, 'orb_beta')
    assert mol.energy == -37.793939631890
    check_orthonormality(mol)


def test_carbon_sc_ae_uncontracted():
    fn_out = get_fn('carbon_sc_ae_uncontracted.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [6]).all()
    assert (mol.pseudo_numbers == [6]).all()
    assert (mol.orb_alpha_occs == [1, 1, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]).all()
    assert (mol.orb_alpha_energies == [-10.062206, -0.499716, -0.192580,
                                       -0.192580, -0.192580]).all()
    assert not hasattr(mol, 'orb_beta')
    assert mol.energy == -37.800453482378
    check_orthonormality(mol)


def test_carbon_sc_pp_contracted():
    fn_out = get_fn('carbon_sc_pp_contracted.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [6]).all()
    assert (mol.pseudo_numbers == [4]).all()
    assert (mol.orb_alpha_occs == [1, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]).all()
    assert (mol.orb_alpha_energies == [-0.500732, -0.193138, -0.193138, -0.193138]).all()
    assert not hasattr(mol, 'orb_beta')
    assert mol.energy == -5.350765755382
    check_orthonormality(mol)


def test_carbon_sc_pp_uncontracted():
    fn_out = get_fn('carbon_sc_pp_uncontracted.cp2k.out')
    mol = IOData.from_file(fn_out)
    assert (mol.numbers == [6]).all()
    assert (mol.pseudo_numbers == [4]).all()
    assert (mol.orb_alpha_occs == [1, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]).all()
    assert (mol.orb_alpha_energies == [-0.500238, -0.192365, -0.192365, -0.192365]).all()
    assert not hasattr(mol, 'orb_beta')
    assert mol.energy == -5.352864672201
    check_orthonormality(mol)


def test_errors():
    fn_test = get_fn('carbon_sc_pp_uncontracted.cp2k.out')
    with truncated_file('io.test.test_cp2k.test_errors', fn_test, 0, 0) as fn:
        with assert_raises(IOError):
            IOData.from_file(fn)
    with truncated_file('io.test.test_cp2k.test_errors', fn_test, 107, 10) as fn:
        with assert_raises(IOError):
            IOData.from_file(fn)
    with truncated_file('io.test.test_cp2k.test_errors', fn_test, 357, 10) as fn:
        with assert_raises(IOError):
            IOData.from_file(fn)
    with truncated_file('io.test.test_cp2k.test_errors', fn_test, 405, 10) as fn:
        with assert_raises(IOError):
            IOData.from_file(fn)
    fn_test = get_fn('carbon_gs_pp_uncontracted.cp2k.out')
    with truncated_file('io.test.test_cp2k.test_errors', fn_test, 456, 10) as fn:
        with assert_raises(IOError):
            IOData.from_file(fn)
