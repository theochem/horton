# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
# Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
#
# This file is part of Horton.
#
# Horton is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# Horton is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--
#pylint: skip-file


from horton import *
from horton.meanfield.test.common import helper_compute


def test_project_identical():
    mol = Molecule.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    exp = mol.lf.create_expansion()
    project_orbitals_mgs(mol.obasis, mol.obasis, mol.exp_alpha, exp)
    assert (exp.energies == 0.0).all()
    assert (exp.occupations == mol.exp_alpha.occupations).all()
    assert abs(exp.coeffs[:,:-2] - mol.exp_alpha.coeffs[:,:-2]).max() < 1e-9
    assert (exp.coeffs[:,-2:] == 0.0).all()


def test_project_larger():
    # Load STO3G system and keep essential results
    mol = Molecule.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
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
    for i0 in xrange(5):
        for i1 in xrange(i0+1):
            dot = olp.inner(exp1.coeffs[:,i0], exp1.coeffs[:,i1])
            if i0 == i1:
                assert abs(dot-1) < 1e-5
            else:
                assert abs(dot) < 1e-5

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


def test_project_smaller():
    # Load 3-21G system and keep essential results
    mol = Molecule.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
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
    for exp, nocc in (exp1_alpha, 2), (exp1_beta, 1):
        for i0 in xrange(nocc):
            for i1 in xrange(i0+1):
                dot = olp.inner(exp.coeffs[:,i0], exp.coeffs[:,i1])
                if i0 == i1:
                    assert abs(dot-1) < 1e-5
                else:
                    assert abs(dot) < 1e-5

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


def test_same_size():
    # Create initial system
    mol = Molecule.from_file(context.get_fn('test/water.xyz'))
    obasis0 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    lf = DenseLinalgFactory(obasis0.nbasis)
    exp0 = lf.create_expansion()


    olp = obasis0.compute_overlap(lf)
    kin = obasis0.compute_kinetic(lf)
    na = obasis0.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    er = obasis0.compute_electron_repulsion(lf)
    guess_core_hamiltonian(olp, kin, na, exp0)

    # Change geometry
    mol.coordinates[1,2] += 0.1
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    exp1 = lf.create_expansion()

    # Project from one to other:
    project_orbitals_mgs(obasis0, obasis1, exp0, exp1)
    assert (exp1.energies == 0.0).all()
    assert (exp1.occupations == exp0.occupations).all()
    assert abs(exp1.coeffs[:,:5] - exp0.coeffs[:,:5]).max() > 1e-3 # something should change
    assert (exp1.coeffs[:,5:] == 0.0).all()
