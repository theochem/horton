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
from nose.tools import assert_raises

from horton.meanfield.test.common import load_olp, load_kin, load_na, load_er, load_nn, load_mdata, \
    load_er_chol
from .. import Orbitals, guess_core_hamiltonian, RTwoIndexTerm, RDirectTerm, RExchangeTerm, REffHam, \
    AufbauOccModel, CDIISSCFSolver
from ..indextransform import _parse_four_index_transform_orbs, four_index_transform, \
    four_index_transform_cholesky, split_core_active, split_core_active_cholesky


def test_parse_index_transform_orbs():
    assert _parse_four_index_transform_orbs(0, 1, 2, 3) == (0, 1, 2, 3)
    assert _parse_four_index_transform_orbs(0, None, None, None) == (0, 0, 0, 0)
    assert _parse_four_index_transform_orbs(0, 1, None, None) == (0, 1, 0, 1)
    assert _parse_four_index_transform_orbs(0, None, 2, None) == (0, 0, 2, 2)
    with assert_raises(TypeError):
        _parse_four_index_transform_orbs(0, None, None, 2)


def test_exceptions():
    with assert_raises(ValueError):
        four_index_transform(np.zeros((4, 4, 4, 4)), 0, method='foo')
    with assert_raises(ValueError):
        four_index_transform_cholesky(np.zeros((2, 4, 4)), 0, method='foo')
    with assert_raises(ValueError):
        split_core_active(np.zeros((5, 5)), None, 0.0, None, -1, 3)
    with assert_raises(ValueError):
        split_core_active(np.zeros((5, 5)), None, 0.0, None, 3, 0)
    with assert_raises(ValueError):
        split_core_active(np.zeros((5, 5)), None, 0.0, None, 3, 7)
    with assert_raises(ValueError):
        split_core_active_cholesky(np.zeros((5, 5)), None, 0.0, None, -1, 3)
    with assert_raises(ValueError):
        split_core_active_cholesky(np.zeros((5, 5)), None, 0.0, None, 3, 0)
    with assert_raises(ValueError):
        split_core_active_cholesky(np.zeros((5, 5)), None, 0.0, None, 3, 7)


def helper_hf(olp, ecore, one, two, nocc):
    # Initial guess
    orb_alpha = Orbitals(olp.shape[0])
    guess_core_hamiltonian(olp, one, orb_alpha)

    # Construct the restricted HF effective Hamiltonian
    external = {'core': ecore}
    terms = [
        RTwoIndexTerm(one, 'one'),
        RDirectTerm(two, 'hartree'),
        RExchangeTerm(two, 'x_hf'),
    ]
    ham = REffHam(terms, external)

    occ_model = AufbauOccModel(nocc)

    # Converge WFN with CDIIS SCF
    occ_model.assign(orb_alpha)
    dm_alpha = orb_alpha.to_dm()
    scf_solver = CDIISSCFSolver(1e-6)
    scf_solver(ham, olp, occ_model, dm_alpha)

    # Reconstruct the orbitals by diagonalizing the final fock matrix
    fock_alpha = np.zeros(olp.shape)
    ham.compute_fock(fock_alpha)
    orb_alpha.from_fock(fock_alpha, olp)

    return ham.cache['energy'], orb_alpha


def prepare_hf(fname):
    # Input structure
    mdata = load_mdata(fname)

    # Compute Gaussian integrals
    olp = load_olp(fname)  # FIXME: Doesn't take into account basis_str
    kin = load_kin(fname)
    na = load_na(fname)
    one = kin + na
    two = load_er(fname)

    enucnuc = load_nn(fname)
    return olp, kin, na, one, two, enucnuc


def check_core_active_blank(fname):
    olp, kin, na, one, two, enucnuc = prepare_hf(fname)
    ncore = 0
    nactive = olp.shape[0]
    one_small, two_small, ecore = \
        split_core_active(one, two, enucnuc, None, ncore, nactive)
    np.testing.assert_allclose(two, two_small)
    np.testing.assert_allclose(enucnuc, ecore)
    np.testing.assert_allclose(one, one_small)


def check_core_active(fname, ncore, nactive):
    # A) Run a simple HF calculation on the given IOData in the given basis
    olp, kin, na, one, two, enucnuc = prepare_hf(fname)

    mdata = load_mdata(fname)
    # Decide how to occupy the orbitals
    assert mdata['numbers'].sum() % 2 == 0
    nocc = mdata['numbers'].sum() / 2
    assert ncore + nactive > nocc

    enucnuc = load_nn(fname)
    energy1, orb_alpha1 = helper_hf(olp, enucnuc, one, two, nocc)

    # B1) Get integrals for the active space, using tensordot transformation
    one_small, two_small, ecore = split_core_active(
        one, two, enucnuc, orb_alpha1, ncore, nactive)
    # C1) Verify the RHF energy using the active space integrals
    energy2, orb_alpha2 = helper_hf(
        np.identity(len(one_small)), ecore, one_small, two_small, nocc - ncore)
    np.testing.assert_almost_equal(energy1, energy2)

    # B2) Get integrals for the active space, using einsum transformation
    one_small, two_small, ecore = split_core_active(
        one, two, enucnuc, orb_alpha1, ncore, nactive, indextrans='einsum')
    # C2) Verify the RHF energy using the active space integrals
    energy2, orb_alpha2 = helper_hf(
        np.identity(len(one_small)), ecore, one_small, two_small, nocc - ncore)
    np.testing.assert_almost_equal(energy1, energy2)


def test_core_active_neon():
    fname = "neon_6_31gd_xyz"
    check_core_active_blank(fname)
    check_core_active(fname, 2, 6)


def test_core_active_water():
    fname = 'water_6_31gd_xyz'
    check_core_active_blank(fname)
    check_core_active(fname, 1, 6)


# TODO: Move to higher level test
# def test_core_active_h2_azirine():
#     fname = 'h2_azirine_3_21g_xyz'
#     check_core_active_blank(fname)
#     check_core_active(fname, 3, 15)


def helper_hf_cholesky(olp, ecore, one, two_vecs, nocc):
    # Initial guess
    orb_alpha = Orbitals(olp.shape[0])
    guess_core_hamiltonian(olp, one, orb_alpha)

    # Construct the restricted HF effective Hamiltonian
    external = {'core': ecore}
    terms = [
        RTwoIndexTerm(one, 'one'),
        RDirectTerm(two_vecs, 'hartree'),
        RExchangeTerm(two_vecs, 'x_hf'),
    ]
    ham = REffHam(terms, external)

    occ_model = AufbauOccModel(nocc)

    # Converge WFN with CDIIS SCF
    occ_model.assign(orb_alpha)
    dm_alpha = orb_alpha.to_dm()
    scf_solver = CDIISSCFSolver(1e-6)
    scf_solver(ham, olp, occ_model, dm_alpha)

    # Reconstruct the orbitals by diagonalizing the final fock matrix
    fock_alpha = np.zeros(olp.shape)
    ham.compute_fock(fock_alpha)
    orb_alpha.from_fock(fock_alpha, olp)

    return ham.cache['energy'], orb_alpha


def prepare_hf_cholesky(fname):
    # Input structure
    mdata = load_mdata(fname)

    # Compute Gaussian integrals
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    one = kin + na
    two_vecs = load_er_chol(fname)

    enucnuc = load_nn(fname)
    return olp, kin, na, one, two_vecs, enucnuc


def check_core_active_cholesky_blank(fname):
    olp, kin, na, one, two_vecs, enucnuc = prepare_hf_cholesky(fname)
    ncore = 0
    nactive = olp.shape[0]
    one_small, two_vecs_small, ecore = \
        split_core_active_cholesky(one, two_vecs, enucnuc, None, ncore, nactive)
    np.testing.assert_allclose(two_vecs, two_vecs_small)
    np.testing.assert_allclose(enucnuc, ecore)
    np.testing.assert_allclose(one, one_small)


def check_core_active_cholesky(fname, ncore, nactive):
    # A) Run a simple HF calculation on the given IOData in the given basis
    olp, kin, na, one, two_vecs, enucnuc = prepare_hf_cholesky(fname)

    mdata = load_mdata(fname)
    # Decide how to occupy the orbitals
    assert mdata['numbers'].sum() % 2 == 0
    nocc = mdata['numbers'].sum() / 2
    assert ncore + nactive > nocc

    energy1, orb_alpha1 = helper_hf_cholesky(olp, enucnuc, one, two_vecs, nocc)

    # B1) Get integrals for the active space, using tensordot transformation
    one_small, two_vecs_small, ecore = split_core_active_cholesky(
        one, two_vecs, enucnuc, orb_alpha1, ncore, nactive)
    # C1) Verify the RHF energy using the active space integrals
    energy2, orb_alpha2 = helper_hf_cholesky(
        np.identity(len(one_small)), ecore, one_small, two_vecs_small, nocc - ncore)
    np.testing.assert_almost_equal(energy1, energy2)

    # B2) Get integrals for the active space, using einsum transformation
    one_small, two_vecs_small, ecore = split_core_active_cholesky(
        one, two_vecs, enucnuc, orb_alpha1, ncore, nactive, indextrans='einsum')
    # C2) Verify the RHF energy using the active space integrals
    energy2, orb_alpha2 = helper_hf_cholesky(
        np.identity(len(one_small)), ecore, one_small, two_vecs_small, nocc - ncore)
    np.testing.assert_almost_equal(energy1, energy2)


def test_core_active_neon_cholesky():
    fname = "neon_6_31gd_xyz"
    check_core_active_cholesky_blank(fname)
    check_core_active_cholesky(fname, 2, 6)


def test_core_active_water_cholesky():
    fname = 'water_6_31gd_xyz'
    check_core_active_cholesky_blank(fname)
    check_core_active_cholesky(fname, 1, 6)

# TODO: Move to higher level test
# def test_core_active_h2_azirine_cholesky():
#     fname = 'h2_azirine_3_21g_xyz'
#     check_core_active_cholesky_blank(fname)
#     check_core_active_cholesky(fname, 3, 15)
