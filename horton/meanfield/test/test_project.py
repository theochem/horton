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


from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import helper_compute


def test_project_msg_identical():
    mol = IOData.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    exp = mol.lf.create_expansion()
    project_orbitals_mgs(mol.obasis, mol.obasis, mol.exp_alpha, exp)
    assert (exp.energies == 0.0).all()
    assert (exp.occupations == mol.exp_alpha.occupations).all()
    assert abs(exp.coeffs[:,:-2] - mol.exp_alpha.coeffs[:,:-2]).max() < 1e-9
    assert (exp.coeffs[:,-2:] == 0.0).all()


def test_project_ortho_basis_identical():
    mol = IOData.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    exp = mol.lf.create_expansion()
    project_orbitals_ortho(mol.obasis, mol.obasis, mol.exp_alpha, exp)
    assert (exp.energies == 0.0).all()
    assert (exp.occupations == mol.exp_alpha.occupations).all()
    assert abs(exp.coeffs - mol.exp_alpha.coeffs).max() < 1e-9


def test_project_ortho_olp_identical():
    mol = IOData.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    olp = mol.lf.create_two_index()
    for i in xrange(olp.nbasis):
        olp.set_element(i, i, 1.0)
    exp = mol.lf.create_expansion()
    project_orbitals_ortho(mol.obasis, mol.obasis, mol.exp_alpha, exp)
    assert (exp.energies == 0.0).all()
    assert (exp.occupations == mol.exp_alpha.occupations).all()
    assert abs(exp.coeffs - mol.exp_alpha.coeffs).max() < 1e-9


def test_project_msg_larger():
    # Load STO3G system and keep essential results
    mol = IOData.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    obasis0 = mol.obasis
    exp0 = mol.exp_alpha

    # Upgrade the basis to 3-21G and project
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, '3-21G')
    lf1 = DenseLinalgFactory(obasis1.nbasis)
    exp1 = lf1.create_expansion()
    project_orbitals_mgs(obasis0, obasis1, exp0, exp1)
    assert (exp1.energies == 0.0).all()
    assert exp0.occupations.sum() == exp1.occupations.sum()
    assert (exp1.coeffs[:,5:] == 0.0).all()

    # Check the normalization of the projected orbitals
    olp = obasis1.compute_overlap(lf1)
    exp1.check_orthonormality(olp)

    # Setup HF hamiltonian and compute energy
    kin = obasis1.compute_kinetic(lf1)
    na = obasis1.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf1)
    er = obasis1.compute_electron_repulsion(lf1)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)

    # Compute energy after projection
    energy1 = helper_compute(ham, lf1, exp1)[0]

    # Optimize wfn
    scf_solver = PlainSCFSolver(1e-6)
    occ_model = AufbauOccModel(5)
    scf_solver(ham, lf1, olp, occ_model, exp1)
    energy2 = ham.cache['energy']
    assert energy2 < energy1 # the energy should decrease after scf convergence

    # Construct a core initial guess
    guess_core_hamiltonian(olp, kin, na, exp1)
    energy3 = helper_compute(ham, lf1, exp1)[0]
    assert energy3 > energy1 # the projected guess should be better than the core guess


def test_project_msg_smaller():
    # Load 3-21G system and keep essential results
    mol = IOData.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
    obasis0 = mol.obasis
    exp0_alpha = mol.exp_alpha
    exp0_beta = mol.exp_beta

    # Downgrade the basis to sto-3g and project
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    lf1 = DenseLinalgFactory(obasis1.nbasis)
    exp1_alpha = lf1.create_expansion()
    exp1_beta = lf1.create_expansion()
    project_orbitals_mgs(obasis0, obasis1, exp0_alpha, exp1_alpha)
    project_orbitals_mgs(obasis0, obasis1, exp0_beta, exp1_beta)
    assert (exp1_alpha.energies == 0.0).all()
    assert (exp1_beta.energies == 0.0).all()
    assert exp1_alpha.occupations.sum() == 2
    assert exp1_beta.occupations.sum() == 1
    assert (exp1_alpha.coeffs[:,2:] == 0.0).all()
    assert (exp1_beta.coeffs[:,1:] == 0.0).all()

    # Check the normalization of the projected orbitals
    olp = obasis1.compute_overlap(lf1)
    exp1_alpha.check_orthonormality(olp)
    exp1_beta.check_orthonormality(olp)

    # Setup HF hamiltonian and compute energy
    kin = obasis1.compute_kinetic(lf1)
    na = obasis1.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf1)
    er = obasis1.compute_electron_repulsion(lf1)
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham = UEffHam(terms)

    # Compute energy before SCF
    energy1 = helper_compute(ham, lf1, exp1_alpha, exp1_beta)[0]

    scf_solver = PlainSCFSolver(1e-6)
    occ_model = AufbauOccModel(2, 1)
    scf_solver(ham, lf1, olp, occ_model, exp1_alpha, exp1_beta)
    energy2 = ham.cache['energy']
    assert energy2 < energy1 # the energy should decrease after scf convergence


def get_basis_pair_geometry():
    '''Prepare two basis sets that only differ in geometry'''
    # Create initial system
    mol = IOData.from_file(context.get_fn('test/water.xyz'))
    obasis0 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    lf = DenseLinalgFactory(obasis0.nbasis)
    exp0 = lf.create_expansion()

    # Occupy all orbitals such that orthogonality is well tested
    exp0.occupations[:] = 1.0

    # core-hamiltonian guess
    olp = obasis0.compute_overlap(lf)
    kin = obasis0.compute_kinetic(lf)
    na = obasis0.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    er = obasis0.compute_electron_repulsion(lf)
    guess_core_hamiltonian(olp, kin, na, exp0)

    # Internal consistency check
    exp0.check_orthonormality(obasis0.compute_overlap(lf))

    # Change geometry
    mol.coordinates[1,2] += 0.5
    mol.coordinates[0,1] -= 1.5
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    exp1 = lf.create_expansion()

    return obasis0, obasis1, exp0, exp1, lf


def test_project_msg_geometry():
    obasis0, obasis1, exp0, exp1, lf = get_basis_pair_geometry()

    # Project from one to other:
    project_orbitals_mgs(obasis0, obasis1, exp0, exp1)

    # Basic checks
    assert (exp1.energies == 0.0).all()
    assert (exp1.occupations == exp0.occupations).all()
    assert abs(exp1.coeffs[:,:5] - exp0.coeffs[:,:5]).max() > 1e-3 # something should change

    # Check orthonormality
    exp1.check_orthonormality(obasis1.compute_overlap(lf))


def test_project_ortho_basis_geometry():
    obasis0, obasis1, exp0, exp1, lf = get_basis_pair_geometry()

    # Project from one to other:
    project_orbitals_ortho(obasis0, obasis1, exp0, exp1)

    # Basic checks
    assert (exp1.energies == 0.0).all()
    assert (exp1.occupations == exp0.occupations).all()
    assert abs(exp1.coeffs[:,:5] - exp0.coeffs[:,:5]).max() > 1e-3 # something should change

    # Check orthonormality
    exp1.check_orthonormality(obasis1.compute_overlap(lf))


def test_project_ortho_olp_geometry():
    obasis0, obasis1, exp0, exp1, lf = get_basis_pair_geometry()

    # Project from one to other:
    olp0 = obasis0.compute_overlap(lf)
    olp1 = obasis1.compute_overlap(lf)
    project_orbitals_ortho(olp0, olp1, exp0, exp1)

    # Basic checks
    assert (exp1.energies == 0.0).all()
    assert (exp1.occupations == exp0.occupations).all()
    assert abs(exp1.coeffs[:,:5] - exp0.coeffs[:,:5]).max() > 1e-3 # something should change

    # Check orthonormality
    exp1.check_orthonormality(obasis1.compute_overlap(lf))
