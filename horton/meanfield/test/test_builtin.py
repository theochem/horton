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

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_becke_hartree_n2_hfs_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False, mode='keep')

    er = mol.obasis.compute_electron_repulsion()
    ham1 = REffHam([RDirectTerm(er, 'hartree')])
    ham2 = REffHam([RGridGroup(mol.obasis, grid, [RBeckeHartree(8)])])

    dm_alpha = mol.orb_alpha.to_dm()
    ham1.reset(dm_alpha)
    ham2.reset(dm_alpha)
    energy1 = ham1.compute_energy()
    energy2 = ham2.compute_energy()
    assert abs(energy1 - energy2) < 1e-3

    op1 = np.zeros(dm_alpha.shape)
    op2 = np.zeros(dm_alpha.shape)
    ham1.compute_fock(op1)
    ham2.compute_fock(op2)
    np.testing.assert_allclose(op1, op2, atol=1e-3)


def test_becke_hartree_h3_hfs_321g():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False, mode='keep')

    er = mol.obasis.compute_electron_repulsion()
    ham1 = UEffHam([UDirectTerm(er, 'hartree')])
    ham2 = UEffHam([UGridGroup(mol.obasis, grid, [UBeckeHartree(8)])])

    dm_alpha = mol.orb_alpha.to_dm()
    dm_beta = mol.orb_beta.to_dm()
    ham1.reset(dm_alpha, dm_beta)
    ham2.reset(dm_alpha, dm_beta)
    energy1 = ham1.compute_energy()
    energy2 = ham2.compute_energy()
    assert abs(energy1 - energy2) < 1e-3

    fock_alpha1 = np.zeros(dm_alpha.shape)
    fock_beta1 = np.zeros(dm_beta.shape)
    fock_alpha2 = np.zeros(dm_alpha.shape)
    fock_beta2 = np.zeros(dm_beta.shape)
    ham1.compute_fock(fock_alpha1, fock_beta1)
    ham2.compute_fock(fock_alpha2, fock_beta2)
    np.testing.assert_allclose(fock_alpha1, fock_alpha2, atol=1e-3)
    np.testing.assert_allclose(fock_beta1, fock_beta2, atol=1e-3)
