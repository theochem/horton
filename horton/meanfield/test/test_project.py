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
from .common import helper_compute


def test_project_msg_identical():
    mol = IOData.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    orb = Orbitals(mol.obasis.nbasis)
    project_orbitals_mgs(mol.obasis, mol.obasis, mol.orb_alpha, orb)
    assert (orb.energies == 0.0).all()
    assert (orb.occupations == mol.orb_alpha.occupations).all()
    assert abs(orb.coeffs[:,:-2] - mol.orb_alpha.coeffs[:,:-2]).max() < 1e-9
    assert (orb.coeffs[:,-2:] == 0.0).all()


def test_project_ortho_basis_identical():
    mol = IOData.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    orb = Orbitals(mol.obasis.nbasis)
    project_orbitals_ortho(mol.obasis, mol.obasis, mol.orb_alpha, orb)
    assert (orb.energies == 0.0).all()
    assert (orb.occupations == mol.orb_alpha.occupations).all()
    assert abs(orb.coeffs - mol.orb_alpha.coeffs).max() < 1e-9


def test_project_ortho_olp_identical():
    mol = IOData.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    olp = np.identity(mol.obasis.nbasis)
    orb = Orbitals(mol.obasis.nbasis)
    project_orbitals_ortho(mol.obasis, mol.obasis, mol.orb_alpha, orb)
    assert (orb.energies == 0.0).all()
    assert (orb.occupations == mol.orb_alpha.occupations).all()
    assert abs(orb.coeffs - mol.orb_alpha.coeffs).max() < 1e-9


def test_project_msg_larger():
    # Load STO3G system and keep essential results
    mol = IOData.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    obasis0 = mol.obasis
    orb0 = mol.orb_alpha

    # Upgrade the basis to 3-21G and project
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, '3-21G')
    orb1 = Orbitals(obasis1.nbasis)
    project_orbitals_mgs(obasis0, obasis1, orb0, orb1)
    assert (orb1.energies == 0.0).all()
    assert orb0.occupations.sum() == orb1.occupations.sum()
    assert (orb1.coeffs[:,5:] == 0.0).all()

    # Check the normalization of the projected orbitals
    olp = obasis1.compute_overlap()
    orb1.check_orthonormality(olp)

    # Setup HF hamiltonian and compute energy
    kin = obasis1.compute_kinetic()
    na = obasis1.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    er = obasis1.compute_electron_repulsion()
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)

    # Compute energy after projection
    energy1 = helper_compute(ham, orb1)[0]

    # Optimize wfn
    scf_solver = PlainSCFSolver(1e-6)
    occ_model = AufbauOccModel(5)
    scf_solver(ham, olp, occ_model, orb1)
    energy2 = ham.cache['energy']
    assert energy2 < energy1 # the energy should decrease after scf convergence

    # Construct a core initial guess
    guess_core_hamiltonian(olp, kin+na, orb1)
    energy3 = helper_compute(ham, orb1)[0]
    assert energy3 > energy1 # the projected guess should be better than the core guess


def test_project_msg_smaller():
    # Load 3-21G system and keep essential results
    mol = IOData.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
    obasis0 = mol.obasis
    orb0_alpha = mol.orb_alpha
    orb0_beta = mol.orb_beta

    # Downgrade the basis to sto-3g and project
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    orb1_alpha = Orbitals(obasis1.nbasis)
    orb1_beta = Orbitals(obasis1.nbasis)
    project_orbitals_mgs(obasis0, obasis1, orb0_alpha, orb1_alpha)
    project_orbitals_mgs(obasis0, obasis1, orb0_beta, orb1_beta)
    assert (orb1_alpha.energies == 0.0).all()
    assert (orb1_beta.energies == 0.0).all()
    assert orb1_alpha.occupations.sum() == 2
    assert orb1_beta.occupations.sum() == 1
    assert (orb1_alpha.coeffs[:,2:] == 0.0).all()
    assert (orb1_beta.coeffs[:,1:] == 0.0).all()

    # Check the normalization of the projected orbitals
    olp = obasis1.compute_overlap()
    orb1_alpha.check_orthonormality(olp)
    orb1_beta.check_orthonormality(olp)

    # Setup HF hamiltonian and compute energy
    kin = obasis1.compute_kinetic()
    na = obasis1.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    er = obasis1.compute_electron_repulsion()
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham = UEffHam(terms)

    # Compute energy before SCF
    energy1 = helper_compute(ham, orb1_alpha, orb1_beta)[0]

    scf_solver = PlainSCFSolver(1e-6)
    occ_model = AufbauOccModel(2, 1)
    scf_solver(ham, olp, occ_model, orb1_alpha, orb1_beta)
    energy2 = ham.cache['energy']
    assert energy2 < energy1 # the energy should decrease after scf convergence


def get_basis_pair_geometry():
    '''Prepare two basis sets that only differ in geometry'''
    # Create initial system
    mol = IOData.from_file(context.get_fn('test/water.xyz'))
    obasis0 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    orb0 = Orbitals(obasis0.nbasis)

    # Occupy all orbitals such that orthogonality is well tested
    orb0.occupations[:] = 1.0

    # core-hamiltonian guess
    olp = obasis0.compute_overlap()
    kin = obasis0.compute_kinetic()
    na = obasis0.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    er = obasis0.compute_electron_repulsion()
    guess_core_hamiltonian(olp, kin+na, orb0)

    # Internal consistency check
    orb0.check_orthonormality(obasis0.compute_overlap())

    # Change geometry
    mol.coordinates[1,2] += 0.5
    mol.coordinates[0,1] -= 1.5
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    orb1 = Orbitals(obasis1.nbasis)

    return obasis0, obasis1, orb0, orb1


def test_project_msg_geometry():
    obasis0, obasis1, orb0, orb1= get_basis_pair_geometry()

    # Project from one to other:
    project_orbitals_mgs(obasis0, obasis1, orb0, orb1)

    # Basic checks
    assert (orb1.energies == 0.0).all()
    assert (orb1.occupations == orb0.occupations).all()
    assert abs(orb1.coeffs[:,:5] - orb0.coeffs[:,:5]).max() > 1e-3 # something should change

    # Check orthonormality
    orb1.check_orthonormality(obasis1.compute_overlap())


def test_project_ortho_basis_geometry():
    obasis0, obasis1, orb0, orb1 = get_basis_pair_geometry()

    # Project from one to other:
    project_orbitals_ortho(obasis0, obasis1, orb0, orb1)

    # Basic checks
    assert (orb1.energies == 0.0).all()
    assert (orb1.occupations == orb0.occupations).all()
    assert abs(orb1.coeffs[:,:5] - orb0.coeffs[:,:5]).max() > 1e-3 # something should change

    # Check orthonormality
    orb1.check_orthonormality(obasis1.compute_overlap())


def test_project_ortho_olp_geometry():
    obasis0, obasis1, orb0, orb1 = get_basis_pair_geometry()

    # Project from one to other:
    olp0 = obasis0.compute_overlap()
    olp1 = obasis1.compute_overlap()
    project_orbitals_ortho(olp0, olp1, orb0, orb1)

    # Basic checks
    assert (orb1.energies == 0.0).all()
    assert (orb1.occupations == orb0.occupations).all()
    assert abs(orb1.coeffs[:,:5] - orb0.coeffs[:,:5]).max() > 1e-3 # something should change

    # Check orthonormality
    orb1.check_orthonormality(obasis1.compute_overlap())
