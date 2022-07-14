# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2022 The HORTON Development Team
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

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_occ_aufbau_cs():
    orb = Orbitals(10)

    # integer
    occ_model = AufbauOccModel(3)
    occ_model.assign(orb)
    assert abs(orb.occupations[:4] - [1.0, 1.0, 1.0, 0.0]).max() < 1e-10

    # fractional 1
    occ_model = AufbauOccModel(2.9)
    occ_model.assign(orb)
    assert abs(orb.occupations[:4] - [1.0, 1.0, 0.9, 0.0]).max() < 1e-10

    # fractional 2
    occ_model = AufbauOccModel(2.8)
    occ_model.assign(orb)
    assert abs(orb.occupations[:4] - [1.0, 1.0, 0.8, 0.0]).max() < 1e-10


def test_occ_aufbau_os():
    orb_alpha = Orbitals(10)
    orb_beta = Orbitals(10)

    # integer
    occ_model = AufbauOccModel(3, 4)
    occ_model.assign(orb_alpha, orb_beta)
    assert abs(orb_alpha.occupations[:5] - [1.0, 1.0, 1.0, 0.0, 0.0]).max() < 1e-10
    assert abs(orb_beta.occupations[:5] - [1.0, 1.0, 1.0, 1.0, 0.0]).max() < 1e-10

    # fractional
    occ_model = AufbauOccModel(2.9, 3.1)
    occ_model.assign(orb_alpha, orb_beta)
    assert abs(orb_alpha.occupations[:5] - [1.0, 1.0, 0.9, 0.0, 0.0]).max() < 1e-10
    assert abs(orb_beta.occupations[:5] - [1.0, 1.0, 1.0, 0.1, 0.0]).max() < 1e-10


def test_fermi_occ_model_cs_helium():
    fn_fchk = context.get_fn('test/helium_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    occ_model = FermiOccModel(1.0)
    occ_model.assign(mol.orb_alpha)
    assert (mol.orb_alpha.occupations == [1.0]).all()


def test_fermi_occ_model_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    for temperature in 300, 3000, 10000, 30000:
        occ_model = FermiOccModel(5.0, temperature=temperature)
        occ_model.assign(mol.orb_alpha)
        occ = mol.orb_alpha.occupations
        assert abs(occ.sum() - 5.0) < 1e-8
        assert (occ[1:] <= occ[:-1]).all()


def test_fermi_occ_model_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = IOData.from_file(fn_fchk)
    for temperature in 300, 3000, 10000, 30000:
        occ_model = FermiOccModel(1.9, 1.1, temperature=temperature)
        occ_model.assign(mol.orb_alpha, mol.orb_beta)
        occ_a = mol.orb_alpha.occupations
        assert abs(occ_a.sum() - 1.9) < 1e-8
        assert (occ_a[1:] <= occ_a[:-1]).all()
        occ_b = mol.orb_beta.occupations
        assert abs(occ_b.sum() - 1.1) < 1e-8
        assert (occ_b[1:] <= occ_b[:-1]).all()


def test_fixed_occ_model_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = IOData.from_file(fn_fchk)
    occs_alpha = np.array([2.0, 0.0, 0.5])
    occs_beta = np.array([0.0, 0.5, 0.0, 0.0])
    occ_model = FixedOccModel(occs_alpha, occs_beta)
    mol.orb_alpha.occupations[:] = 0.2
    occ_model.assign(mol.orb_alpha, mol.orb_beta)
    assert (mol.orb_alpha.occupations[:len(occs_alpha)] == occs_alpha).all()
    assert (mol.orb_alpha.occupations[len(occs_alpha):] == 0.0).all()
    assert (mol.orb_beta.occupations[:len(occs_beta)] == occs_beta).all()
    assert (mol.orb_beta.occupations[len(occs_beta):] == 0.0).all()
