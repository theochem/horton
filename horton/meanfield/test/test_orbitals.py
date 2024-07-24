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
import h5py as h5
import pytest

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.test.common import numpy_seed


#
# Utility functions
#


def get_forth_back(n):
    """Return matching pair of forth and back permutation.

    Parameters
    ----------
    n : int
        The length of the permutation.

    It is guaranteed that the identity permutation is never returned.
    """
    while True:
        forth = np.random.uniform(0, 1, 5).argsort()
        if (forth != np.arange(5)).all():
            break
    back = np.zeros(5, int)
    for i, j in enumerate(forth):
        back[j] = i
    return forth, back


def get_signs(n):
    """Return an array with signs (all elements are just +1 or -1).

    Parameters
    ----------
    n : int
        The length of the permutation.

    It is guaranteed that not all signs are positive.
    """
    while True:
        signs = np.random.randint(0, 1, n) * 2 - 1
        if (signs < 0).any():
            return signs


def get_random_orbitals(nbasis):
    """Return a random expansion and an identity overlap matrix."""
    orb = Orbitals(nbasis)
    a = np.random.normal(0, 1, (nbasis, nbasis))
    a = a + a.T
    evals, evecs = np.linalg.eigh(a)
    orb.coeffs[:] = evecs
    orb.occupations[:nbasis // 2] = 1.0
    orb.energies[:] = np.random.uniform(-1, 1, nbasis)
    orb.energies.sort()
    olp = np.identity(nbasis)
    return orb, olp


#
# Actual tests
#

def test_orbitals_hdf5():
    for args in (4,), (6, 3):
        a = Orbitals(*args)
        a.randomize()
        with h5.File('horton.meanfield.test.test_orbitals.test_orbitals_hdf5', "w",
                     driver='core', backing_store=False) as f:
            a.to_hdf5(f)
            b = Orbitals.from_hdf5(f)
            assert a == b
            f.attrs['class'] = 'bork'
            with pytest.raises(TypeError):
                Orbitals.from_hdf5(f)


def test_orbitals_copy_new_randomize_clear_assign():
    for args in (4,), (6, 3):
        a = Orbitals(*args)
        b = a.copy()
        b.randomize()
        assert a != b
        c = b.copy()
        assert b == c
        assert not (b is c)
        d = Orbitals(*args)
        assert a == d
        b.clear()
        assert a == b
        b.assign(c)
        assert b == c


def test_orbitals_copy():
    orb1 = Orbitals(3, 2)
    orb1._coeffs[:] = np.random.uniform(0, 1, (3, 2))
    orb1._energies[:] = np.random.uniform(0, 1, 2)
    orb1._occupations[:] = np.random.uniform(0, 1, 2)
    orb2 = orb1.copy()
    assert (orb1._coeffs == orb2._coeffs).all()
    assert (orb1._energies == orb2._energies).all()
    assert (orb1._occupations == orb2._occupations).all()


def test_orbitals_permute_basis():
    for i in range(10):
        forth, back = get_forth_back(5)
        a = Orbitals(5)
        a.randomize()
        b = a.copy()
        b.permute_basis(forth)
        assert a != b
        b.permute_basis(back)
        assert a == b


def test_orbitals_permute_orbitals():
    for i in range(10):
        forth, back = get_forth_back(5)
        a = Orbitals(5)
        a.randomize()
        b = a.copy()
        b.permute_orbitals(forth)
        assert a != b
        b.permute_orbitals(back)
        assert a == b


def test_orbitals_change_basis_signs():
    for i in range(10):
        signs = get_signs(5)
        a = Orbitals(5)
        a.randomize()
        b = a.copy()
        b.change_basis_signs(signs)
        assert a != b
        b.change_basis_signs(signs)
        assert a == b


def test_orbitals_check_normalization():
    orb, olp = get_random_orbitals(5)
    orb.check_normalization(olp)


def test_orbitals_check_orthonormality():
    orb, olp = get_random_orbitals(5)
    orb.occupations[0] = 0.0
    orb.occupations[-1] = 1.0
    orb.check_orthonormality(olp)


def test_orbitals_error_eigen():
    with numpy_seed(1):
        orb = Orbitals(5)
        a = np.random.normal(0, 1, (5, 5))
        fock = a + a.T
        evals, evecs = np.linalg.eigh(fock)
        orb.coeffs[:] = evecs
        orb.energies[:] = evals
        olp = np.identity(5)
        assert orb.error_eigen(fock, olp) < 1e-10
        orb.coeffs[:] += np.random.normal(0, 1e-3, (5, 5))
        assert orb.error_eigen(fock, olp) > 1e-10


def test_orbitals_from_fock():
    with numpy_seed(1):
        a = np.random.normal(0, 1, (5, 5))
        fock = a + a.T
        a = np.random.normal(0, 1, (5, 5))
        olp = np.dot(a, a.T)
        orb = Orbitals(5)
        orb.from_fock(fock, olp)
        assert orb.error_eigen(fock, olp) < 1e-5


def test_orbitals_from_fock_and_dm():
    natom = 5

    # Use a simple Huckel-like model to construct degenerate levels
    fock = np.zeros((natom, natom))
    olp = np.zeros((natom, natom))
    for i in range(natom):
        fock[i, i] = 0.6
        fock[i, (i + 1) % natom] = -0.2
        fock[(i + 1) % natom, i] = -0.2
        olp[i, i] = 1.0
        olp[i, (i + 1) % natom] = 0.2
        olp[(i + 1) % natom, i] = 0.2

    # Create orbitals that will be used to construct various density matrices
    orb = Orbitals(natom)
    orb.from_fock(fock, olp)

    # Checks for every case
    def check_case(orb0):
        dm = orb0.to_dm()
        orb1 = Orbitals(natom)
        orb1.from_fock_and_dm(fock, dm, olp)
        np.testing.assert_almost_equal(orb0.occupations, orb1.occupations)
        assert orb1.error_eigen(fock, olp) < 1e-5
        sds = np.dot(olp, np.dot(dm, olp))
        orb1.energies[:] = orb1.occupations
        assert orb1.error_eigen(sds, olp) < 1e-5

    # Case 1: not difficult, i.e. compatible degeneracies
    orb.occupations[:] = [1, 1, 1, 0, 0]
    check_case(orb)

    # Case 2: incompatible degeneracies
    orb.occupations[:] = [2, 2, 1, 0, 0]
    check_case(orb)

    # Case 3: incompatible degeneracies and rotated degenerate orbitals
    orb.occupations[:] = [2, 1, 0, 0, 0]
    for i in range(36):
        orb.rotate_2orbitals(np.pi / 18.0, 1, 2)
        check_case(orb)

    # Case 4: incompatible degeneracies, fractional occupations and rotated
    # degenerate orbitals
    orb.occupations[:] = [1.5, 0.7, 0.3, 0, 0]
    for i in range(36):
        orb.rotate_2orbitals(np.pi / 18.0, 1, 2)
        check_case(orb)


def test_orbitals_naturals():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    overlap = mol.obasis.compute_overlap()
    dm = mol.orb_alpha.to_dm()
    orb = Orbitals(dm.shape[0])
    orb.derive_naturals(dm, overlap)
    assert orb.occupations.min() > -1e-6
    assert orb.occupations.max() < 1 + 1e-6
    orb.check_normalization(overlap)


def test_orbitals_homo_lumo_ch3_hf():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    assert mol.orb_alpha.get_homo_index() == 4
    assert mol.orb_beta.get_homo_index() == 3
    assert mol.orb_alpha.get_lumo_index() == 5
    assert mol.orb_beta.get_lumo_index() == 4
    assert mol.orb_alpha.get_homo_index(1) == 3
    assert mol.orb_beta.get_homo_index(1) == 2
    assert mol.orb_alpha.get_lumo_index(1) == 6
    assert mol.orb_beta.get_lumo_index(1) == 5
    assert abs(mol.orb_alpha.get_homo_energy() - -3.63936540E-01) < 1e-8
    assert abs(mol.orb_alpha.get_homo_energy(1) - -5.37273275E-01) < 1e-8
    assert abs(mol.orb_alpha.get_lumo_energy() - 6.48361367E-01) < 1e-8
    assert abs(mol.orb_beta.get_homo_energy() - -5.18988806E-01) < 1e-8
    assert abs(mol.orb_beta.get_homo_energy(1) - -5.19454722E-01) < 1e-8
    assert abs(mol.orb_beta.get_lumo_energy() - 3.28562907E-01) < 1e-8
    assert abs(mol.orb_alpha.homo_energy - -3.63936540E-01) < 1e-8
    assert abs(mol.orb_alpha.lumo_energy - 6.48361367E-01) < 1e-8
    assert abs(mol.orb_beta.homo_energy - -5.18988806E-01) < 1e-8
    assert abs(mol.orb_beta.lumo_energy - 3.28562907E-01) < 1e-8
    with pytest.raises(ValueError):
        mol.orb_alpha.get_homo_index(-1)
    with pytest.raises(ValueError):
        mol.orb_alpha.get_lumo_index(-1)


def test_orbitals_to_dm1():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = IOData.from_file(fn_fchk)
    dm = mol.orb_alpha.to_dm()
    dm *= 2
    np.testing.assert_almost_equal(dm, mol.get_dm_full())
    np.testing.assert_almost_equal(dm, dm.T)


def test_orbitals_to_dm2():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    dm = mol.orb_alpha.to_dm() + mol.orb_beta.to_dm()
    olp = mol.obasis.compute_overlap()
    np.testing.assert_almost_equal(dm, mol.get_dm_full())
    np.testing.assert_almost_equal(dm, dm.T)


def test_orbitals_to_dm3():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = IOData.from_file(fn_fchk)
    dm = mol.orb_alpha.to_dm(other=mol.orb_beta)
    assert (dm != dm.T).any()


def test_orbitals_rotate_random():
    orb0, olp = get_random_orbitals(5)
    orb0.check_normalization(olp)
    orb1 = orb0.copy()
    orb1.rotate_random()
    orb1.check_normalization(olp)
    dots = np.dot(orb0.coeffs.T, orb1.coeffs)
    assert not np.allclose(dots, np.identity(5))


def test_orbitals_two_index_rotate_2orbitals():
    orb0, olp = get_random_orbitals(4)
    orb0.check_normalization(olp)
    orb1 = orb0.copy()
    orb1.rotate_2orbitals()
    orb1.check_normalization(olp)
    check = np.identity(4, float)
    dots = np.dot(orb0.coeffs.T, orb1.coeffs)
    check = np.identity(4)
    check[1, 1] = 1.0 / np.sqrt(2)
    check[1, 2] = 1.0 / np.sqrt(2)
    check[2, 1] = -1.0 / np.sqrt(2)
    check[2, 2] = 1.0 / np.sqrt(2)
    np.testing.assert_almost_equal(dots, check)


def test_orbitals_swap_orbitals():
    orb0, olp = get_random_orbitals(4)
    orb0.check_normalization(olp)
    orb1 = orb0.copy()
    orb1.swap_orbitals(np.array([[0, 1], [2, 3]]))
    dots = np.dot(orb0.coeffs.T, orb1.coeffs)
    check = np.zeros((4, 4))
    check[0, 1] = 1.0
    check[1, 0] = 1.0
    check[2, 3] = 1.0
    check[3, 2] = 1.0
    np.testing.assert_almost_equal(dots, check)
    with pytest.raises(TypeError):
        orb1.swap_orbitals(np.zeros((3, 3), dtype=int))
