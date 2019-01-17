# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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
    lf = DenseLinalgFactory(10)
    exp = lf.create_expansion()

    # integer
    occ_model = AufbauOccModel(3)
    occ_model.assign(exp)
    assert abs(exp.occupations[:4] - [1.0, 1.0, 1.0, 0.0]).max() < 1e-10

    # fractional 1
    occ_model = AufbauOccModel(2.9)
    occ_model.assign(exp)
    assert abs(exp.occupations[:4] - [1.0, 1.0, 0.9, 0.0]).max() < 1e-10

    # fractional 2
    occ_model = AufbauOccModel(2.8)
    occ_model.assign(exp)
    assert abs(exp.occupations[:4] - [1.0, 1.0, 0.8, 0.0]).max() < 1e-10


def test_occ_aufbau_os():
    lf = DenseLinalgFactory(10)
    exp_alpha = lf.create_expansion()
    exp_beta = lf.create_expansion()

    # integer
    occ_model = AufbauOccModel(3, 4)
    occ_model.assign(exp_alpha, exp_beta)
    assert abs(exp_alpha.occupations[:5] - [1.0, 1.0, 1.0, 0.0, 0.0]).max() < 1e-10
    assert abs(exp_beta.occupations[:5] - [1.0, 1.0, 1.0, 1.0, 0.0]).max() < 1e-10

    # fractional
    occ_model = AufbauOccModel(2.9, 3.1)
    occ_model.assign(exp_alpha, exp_beta)
    assert abs(exp_alpha.occupations[:5] - [1.0, 1.0, 0.9, 0.0, 0.0]).max() < 1e-10
    assert abs(exp_beta.occupations[:5] - [1.0, 1.0, 1.0, 0.1, 0.0]).max() < 1e-10


def test_fermi_occ_model_cs_helium():
    fn_fchk = context.get_fn('test/helium_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    occ_model = FermiOccModel(1.0)
    occ_model.assign(mol.exp_alpha)
    assert (mol.exp_alpha.occupations == [1.0]).all()


def test_fermi_occ_model_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    for temperature in 300, 3000, 10000, 30000:
        occ_model = FermiOccModel(5.0, temperature=temperature)
        occ_model.assign(mol.exp_alpha)
        occ = mol.exp_alpha.occupations
        assert abs(occ.sum() - 5.0) < 1e-8
        assert (occ[1:] <= occ[:-1]).all()


def test_fermi_occ_model_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = IOData.from_file(fn_fchk)
    for temperature in 300, 3000, 10000, 30000:
        occ_model = FermiOccModel(1.9, 1.1, temperature=temperature)
        occ_model.assign(mol.exp_alpha, mol.exp_beta)
        occ_a = mol.exp_alpha.occupations
        assert abs(occ_a.sum() - 1.9) < 1e-8
        assert (occ_a[1:] <= occ_a[:-1]).all()
        occ_b = mol.exp_beta.occupations
        assert abs(occ_b.sum() - 1.1) < 1e-8
        assert (occ_b[1:] <= occ_b[:-1]).all()


def test_fixed_occ_model_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = IOData.from_file(fn_fchk)
    occs_alpha = np.array([2.0, 0.0, 0.5])
    occs_beta = np.array([0.0, 0.5, 0.0, 0.0])
    occ_model = FixedOccModel(occs_alpha, occs_beta)
    mol.exp_alpha.occupations[:] = 0.2
    occ_model.assign(mol.exp_alpha, mol.exp_beta)
    assert (mol.exp_alpha.occupations[:len(occs_alpha)] == occs_alpha).all()
    assert (mol.exp_alpha.occupations[len(occs_alpha):] == 0.0).all()
    assert (mol.exp_beta.occupations[:len(occs_beta)] == occs_beta).all()
    assert (mol.exp_beta.occupations[len(occs_beta):] == 0.0).all()
