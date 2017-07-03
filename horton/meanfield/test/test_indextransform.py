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


def check_core_active(mol, basis_str, ncore, nactive):
    # A) Run a simple HF calculation on the given IOData in the given basis

    # Input structure
    obasis = get_gobasis(mol.coordinates, mol.numbers, basis_str)

    # Compute Gaussian integrals
    olp = obasis.compute_overlap()
    kin = obasis.compute_kinetic()
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    one = kin + na
    two = obasis.compute_electron_repulsion()

    # Decide how to occupy the orbitals
    assert mol.numbers.sum() % 2 == 0
    nocc = mol.numbers.sum()/2
    assert ncore + nactive > nocc

    enucnuc = compute_nucnuc(mol.coordinates, mol.pseudo_numbers)
    energy1, orb_alpha1 = helper_hf(olp, enucnuc, one, two, nocc)

    # B1) Get integrals for the active space, using tensordot transformation
    one_small, two_small, ecore = split_core_active(
        one, two, enucnuc, orb_alpha1, ncore, nactive)
    # C1) Verify the RHF energy using the active space integrals
    energy2, orb_alpha2 = helper_hf(
        np.identity(len(one_small)), ecore, one_small, two_small, nocc-ncore)
    np.testing.assert_almost_equal(energy1, energy2)

    # B2) Get integrals for the active space, using einsum transformation
    one_small, two_small, ecore = split_core_active(
        one, two, enucnuc, orb_alpha1, ncore, nactive, indextrans='einsum')
    # C2) Verify the RHF energy using the active space integrals
    energy2, orb_alpha2 = helper_hf(
        np.identity(len(one_small)), ecore, one_small, two_small, nocc-ncore)
    np.testing.assert_almost_equal(energy1, energy2)


def test_core_active_neon():
    mol = IOData(
        coordinates=np.zeros((1, 3), float),
        numbers=np.array([10], int)
    )
    check_core_active(mol, '6-31+g(d)', 2, 6)


def test_core_active_water():
    mol = IOData.from_file(context.get_fn('test/water.xyz'))
    check_core_active(mol, '6-31+g(d)', 1, 6)


def test_core_active_2h_azirine():
    mol = IOData.from_file(context.get_fn('test/2h-azirine.xyz'))
    check_core_active(mol, '3-21g', 3, 15)


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


def check_core_active_cholesky(mol, basis_str, ncore, nactive):
    # A) Run a simple HF calculation on the given IOData in the given basis

    # Input structure
    obasis = get_gobasis(mol.coordinates, mol.numbers, basis_str)

    # Compute Gaussian integrals
    olp = obasis.compute_overlap()
    kin = obasis.compute_kinetic()
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    one = kin + na
    two_vecs = obasis.compute_electron_repulsion_cholesky()

    # Decide how to occupy the orbitals
    assert mol.numbers.sum() % 2 == 0
    nocc = mol.numbers.sum()/2
    assert ncore + nactive > nocc

    enucnuc = compute_nucnuc(mol.coordinates, mol.pseudo_numbers)
    energy1, orb_alpha1 = helper_hf_cholesky(olp, enucnuc, one, two_vecs, nocc)

    # B1) Get integrals for the active space, using tensordot transformation
    one_small, two_vecs_small, ecore = split_core_active_cholesky(
        one, two_vecs, enucnuc, orb_alpha1, ncore, nactive)
    # C1) Verify the RHF energy using the active space integrals
    energy2, orb_alpha2 = helper_hf_cholesky(
        np.identity(len(one_small)), ecore, one_small, two_vecs_small, nocc-ncore)
    np.testing.assert_almost_equal(energy1, energy2)

    # B2) Get integrals for the active space, using einsum transformation
    one_small, two_vecs_small, ecore = split_core_active_cholesky(
        one, two_vecs, enucnuc, orb_alpha1, ncore, nactive, indextrans='einsum')
    # C2) Verify the RHF energy using the active space integrals
    energy2, orb_alpha2 = helper_hf_cholesky(
        np.identity(len(one_small)), ecore, one_small, two_vecs_small, nocc-ncore)
    np.testing.assert_almost_equal(energy1, energy2)


def test_core_active_neon_cholesky():
    mol = IOData(
        coordinates=np.zeros((1, 3), float),
        numbers=np.array([10], int)
    )
    check_core_active_cholesky(mol, '6-31+g(d)', 2, 6)


def test_core_active_water_cholesky():
    mol = IOData.from_file(context.get_fn('test/water.xyz'))
    check_core_active_cholesky(mol, '6-31+g(d)', 1, 6)


def test_core_active_2h_azirine_cholesky():
    mol = IOData.from_file(context.get_fn('test/2h-azirine.xyz'))
    check_core_active_cholesky(mol, '3-21g', 3, 15)
