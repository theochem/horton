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

def check_core_active(mol, basis_str, ncore, nactive):
    #
    # Run a simple HF calculation on the given IOData in the given basis
    #

    # Input structure
    obasis = get_gobasis(mol.coordinates, mol.numbers, basis_str)
    lf = DenseLinalgFactory(obasis.nbasis)

    # Compute Gaussian integrals
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    er = obasis.compute_electron_repulsion(lf)

    # Initial guess
    exp_alpha = lf.create_expansion()
    guess_core_hamiltonian(olp, kin, na, exp_alpha)

    # Construct the restricted HF effective Hamiltonian
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms, external)

    # Decide how to occupy the orbitals
    assert mol.numbers.sum()%2 == 0
    nocc = mol.numbers.sum()/2
    assert ncore + nactive > nocc
    occ_model = AufbauOccModel(nocc)

    # Converge WFN with plain SCF
    occ_model.assign(exp_alpha)
    dm_alpha = exp_alpha.to_dm()
    scf_solver = CDIISSCFSolver(1e-6)
    scf_solver(ham, lf, olp, occ_model, dm_alpha)

    # Reconstruct the orbitals by diagonalizing the final fock matrix
    fock_alpha = lf.create_two_index()
    ham.compute_fock(fock_alpha)
    exp_alpha.from_fock(fock_alpha, olp)

    #
    # Get integrals for the active space
    #

    # Sum all one-body stuff
    one = kin.copy()
    one.iadd(na)

    # Do the actual work that needs to be tested
    one_small, two_small, ecore = split_core_active(one, er, external['nn'], exp_alpha, ncore, nactive)

    #
    # Verify (parts of) the RHF energy using the active space integrals
    #

    # Get the integrals in the mo basis
    (one_mo,), (two_mo,) = transform_integrals(one, er, 'tensordot', exp_alpha)

    # Check the core energy
    ecore_check = external['nn'] + 2*one_mo.trace(0, ncore, 0, ncore)
    ecore_check += two_mo.slice_to_two('abab->ab', None, 2.0, True, 0, ncore, 0, ncore, 0, ncore, 0, ncore).sum()
    ecore_check += two_mo.slice_to_two('abba->ab', None,-1.0, True, 0, ncore, 0, ncore, 0, ncore, 0, ncore).sum()
    assert abs(ecore - ecore_check) < 1e-10

    # Check the one-body energy of the active space
    nocc_small = nocc - ncore
    e_one_active = 2*one_small.trace(0, nocc_small, 0, nocc_small)
    e_one_active_check = 2*one_mo.trace(ncore, nocc, ncore, nocc)
    e_one_active_check += 2*two_mo.slice_to_two('abab->ab', None, 2.0, True,
                                                0, ncore, ncore, nocc,
                                                0, ncore, ncore, nocc).sum()
    e_one_active_check += 2*two_mo.slice_to_two('abba->ab', None,-1.0, True,
                                                0, ncore, ncore, nocc,
                                                ncore, nocc, 0, ncore).sum()
    assert abs(e_one_active - e_one_active_check) < 1e-10

    # Check the two-body energy of the active space
    e_two_active = \
        two_small.slice_to_two('abab->ab', None, 2.0, True, 0, nocc_small, 0, nocc_small, 0, nocc_small, 0, nocc_small).sum() \
      + two_small.slice_to_two('abba->ab', None,-1.0, True, 0, nocc_small, 0, nocc_small, 0, nocc_small, 0, nocc_small).sum()
    e_two_active_check = \
        two_mo.slice_to_two('abab->ab', None, 2.0, True, ncore, nocc, ncore, nocc, ncore, nocc, ncore, nocc).sum() \
      + two_mo.slice_to_two('abba->ab', None,-1.0, True, ncore, nocc, ncore, nocc, ncore, nocc, ncore, nocc).sum()
    assert abs(e_two_active - e_two_active_check) < 1e-10

    # Check the total RHF energy
    e_rhf = ecore + e_one_active + e_two_active
    assert abs(e_rhf - ham.cache['energy']) < 1e-10


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
