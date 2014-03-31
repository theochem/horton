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


def test_project_identical():
    mol = Molecule.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    exp = mol.lf.create_expansion()
    project_orbitals_mgs_low(mol.obasis, mol.obasis, mol.wfn.exp_alpha, exp)
    assert (exp.energies == 0.0).all()
    assert (exp.occupations == mol.wfn.exp_alpha.occupations).all()
    assert abs(exp.coeffs[:,:-2] - mol.wfn.exp_alpha.coeffs[:,:-2]).max() < 1e-9
    assert (exp.coeffs[:,-2:] == 0.0).all()


def test_project_larger():
    # Load STO3G system and keep essential results
    mol = Molecule.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    obasis0 = mol.obasis
    wfn0 = mol.wfn
    exp0 = wfn0.exp_alpha

    # Upgrade the basis to 3-21G and project
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, '3-21G')
    mol.lf.set_default_nbasis(obasis1.nbasis)
    wfn1 = setup_mean_field_wfn(obasis1.nbasis, mol.pseudo_numbers, mol.lf, restricted=True)
    project_orbitals_mgs(obasis0, obasis1, wfn0, wfn1)
    exp1 = wfn1.exp_alpha
    assert (exp1.energies == 0.0).all()
    assert exp0.occupations.sum() == exp1.occupations.sum()
    assert (exp1.coeffs[:,5:] == 0.0).all()

    # Check the normalization of the projected orbitals
    olp = obasis1.compute_overlap(mol.lf)
    for i0 in xrange(5):
        for i1 in xrange(i0+1):
            dot = olp.dot(wfn1.exp_alpha.coeffs[:,i0], wfn1.exp_alpha.coeffs[:,i1])
            if i0 == i1:
                assert abs(dot-1) < 1e-5
            else:
                assert abs(dot) < 1e-5

    # Setup HF hamiltonian and compute energy
    kin = obasis1.compute_kinetic(mol.lf)
    nai = obasis1.compute_nuclear_attraction(mol.pseudo_numbers, mol.coordinates, mol.lf)
    er = obasis1.compute_electron_repulsion(mol.lf)
    terms = [
        OneBodyTerm(kin, mol.lf, wfn1, 'kin'),
        DirectTerm(er, mol.lf, wfn1),
        ExchangeTerm(er, mol.lf, wfn1),
        OneBodyTerm(nai, mol.lf, wfn1, 'ne'),
    ]
    ham = Hamiltonian(terms)
    energy1 = ham.compute()

    # Optimize wfn
    converge_scf_oda(ham, wfn1, mol.lf, olp)
    energy2 = ham.cache['energy']
    assert energy2 < energy1 # the energy should decrease after scf convergence

    # Construct a core initial guess
    guess_core_hamiltonian(wfn1, olp, kin, nai)
    ham.clear()
    energy3 = ham.compute()
    assert energy3 > energy1 # the projected guess should be better than the core guess


def test_project_smaller():
    # Load 3-21G system and keep essential results
    mol = Molecule.from_file(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))
    obasis0 = mol.obasis
    wfn0 = mol.wfn

    # Downgrade the basis to sto-3g and project
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    mol.lf.set_default_nbasis(obasis1.nbasis)
    wfn1 = setup_mean_field_wfn(obasis1.nbasis, mol.numbers, mol.lf, restricted=False)
    project_orbitals_mgs(obasis0, obasis1, wfn0, wfn1)
    assert (wfn1.exp_alpha.energies == 0.0).all()
    assert (wfn1.exp_beta.energies == 0.0).all()
    assert wfn1.exp_alpha.occupations.sum() == 2
    assert wfn1.exp_beta.occupations.sum() == 1
    assert (wfn1.exp_alpha.coeffs[:,2:] == 0.0).all()
    assert (wfn1.exp_beta.coeffs[:,1:] == 0.0).all()

    # Check the normalization of the projected orbitals
    olp = obasis1.compute_overlap(mol.lf)
    for exp, nocc in (wfn1.exp_alpha, 2), (wfn1.exp_beta, 1):
        for i0 in xrange(nocc):
            for i1 in xrange(i0+1):
                dot = olp.dot(exp.coeffs[:,i0], exp.coeffs[:,i1])
                if i0 == i1:
                    assert abs(dot-1) < 1e-5
                else:
                    assert abs(dot) < 1e-5

    # Setup HF hamiltonian and compute energy
    kin = obasis1.compute_kinetic(mol.lf)
    nai = obasis1.compute_nuclear_attraction(mol.pseudo_numbers, mol.coordinates, mol.lf)
    er = obasis1.compute_electron_repulsion(mol.lf)
    terms = [
        OneBodyTerm(kin, mol.lf, wfn1, 'kin'),
        DirectTerm(er, mol.lf, wfn1),
        ExchangeTerm(er, mol.lf, wfn1),
        OneBodyTerm(nai, mol.lf, wfn1, 'ne'),
    ]
    ham = Hamiltonian(terms)

    energy1 = ham.compute()

    # Optimize wfn
    converge_scf_oda(ham, wfn1, mol.lf, olp)
    energy2 = ham.cache['energy']
    assert energy2 < energy1 # the energy should decrease after scf convergence

    # Construct a core initial guess
    guess_core_hamiltonian(wfn1, olp, kin, nai)
    ham.clear()
    energy3 = ham.compute()
    assert energy3 > energy2 # the core guess should be worse than the converged


def test_same_size():
    # Create initial system
    mol = Molecule.from_file(context.get_fn('test/water.xyz'))
    obasis0 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    lf = DenseLinalgFactory(obasis0.nbasis)
    wfn0 = setup_mean_field_wfn(obasis0.nbasis, mol.pseudo_numbers, lf, restricted=True)
    olp = obasis0.compute_overlap(mol.lf)
    kin = obasis0.compute_kinetic(mol.lf)
    nai = obasis0.compute_nuclear_attraction(mol.pseudo_numbers, mol.coordinates, mol.lf)
    er = obasis0.compute_electron_repulsion(mol.lf)
    guess_core_hamiltonian(wfn0, olp, kin, nai)

    # Change geometry
    mol.coordinates[1,2] += 0.1
    obasis1 = get_gobasis(mol.coordinates, mol.numbers, 'sto-3g')
    wfn1 = setup_mean_field_wfn(obasis1.nbasis, mol.pseudo_numbers, lf, restricted=True)

    # Project from one to other:
    project_orbitals_mgs(obasis0, obasis1, wfn0, wfn1)
    assert (wfn1.exp_alpha.energies == 0.0).all()
    assert (wfn1.exp_alpha.occupations == wfn0.exp_alpha.occupations).all()
    assert abs(wfn1.exp_alpha.coeffs[:,:5] - wfn0.exp_alpha.coeffs[:,:5]).max() > 1e-3 # something should change
    assert (wfn1.exp_alpha.coeffs[:,5:] == 0.0).all()
