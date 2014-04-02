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


import numpy as np
from nose.tools import assert_raises

from horton import *
from horton.meanfield.test.common import check_cubic_cs_wrapper


def test_energy_hydrogen():
    fn_fchk = context.get_fn('test/h_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        UnrestrictedOneBodyTerm(kin, 'kin'),
        UnrestrictedDirectTerm(er, 'hartree'),
        UnrestrictedExchangeTerm(er, 'x_hf'),
        UnrestrictedOneBodyTerm(na, 'ne'),
    ]
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    ham = UnrestrictedEffectiveHamiltonian(terms, external)
    ham.reset(mol.wfn.dm_alpha, mol.wfn.dm_beta)
    ham.compute()
    assert abs(ham.cache['energy'] - -4.665818503844346E-01) < 1e-8


def test_energy_n2_hfs_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RestrictedOneBodyTerm(kin, 'kin'),
        RestrictedDirectTerm(er, 'hartree'),
        RestrictedGridGroup(mol.obasis, grid, [
            RestrictedDiracExchange(),
        ]),
        RestrictedOneBodyTerm(na, 'ne'),
    ]
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    ham = RestrictedEffectiveHamiltonian(terms, external)
    ham.reset(mol.wfn.dm_alpha)
    ham.compute()

    # Compare energies
    assert abs(ham.cache['energy_ne'] - -2.981579553570E+02) < 1e-6
    assert abs(ham.cache['energy_kin'] - 1.061620887711E+02) < 1e-6
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_x_dirac'] - 6.247259253877E+01) < 1e-4
    assert abs(ham.cache['energy'] - -106.205213597) < 1e-4
    assert abs(ham.cache['energy_nn'] - 23.3180604505) < 3e-9


    # Test if the grid potential data can be properly converted into an operator:
    pot = ham.cache.load('pot_x_dirac_alpha')
    ev1 = grid.integrate(pot, ham.cache.load('rho_alpha'))
    op = mol.lf.create_one_body()
    mol.obasis.compute_grid_density_fock(grid.points, grid.weights, pot, op)
    ev2 = op.expectation_value(mol.wfn.dm_alpha)
    assert abs(ev1 - ev2) < 1e-10
    # check symmetry
    op.check_symmetry()

    # When repeating, we should get the same
    ham.compute()
    assert abs(ham.cache['energy'] - -106.205213597) < 1e-4


def test_fock_n2_hfs_sto3g():
    # The fock operator is tested by doing an SCF and checking the converged
    # energies
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'veryfine', random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RestrictedOneBodyTerm(kin, 'kin'),
        RestrictedDirectTerm(er, 'hartree'),
        RestrictedGridGroup(mol.obasis, grid, [
            RestrictedDiracExchange(),
        ]),
        RestrictedOneBodyTerm(na, 'ne'),
    ]
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    ham = RestrictedEffectiveHamiltonian(terms, external)

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 1e-5

    # Converge from scratch
    guess_core_hamiltonian(mol.wfn, olp, kin, na)
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) > 1e-8
    converge_scf(ham, mol.wfn, mol.lf, olp)
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 1e-8

    # test orbital energies
    expected_energies = np.array([
        -1.37107053E+01, -1.37098006E+01, -9.60673085E-01, -3.57928483E-01,
        -3.16017655E-01, -3.16017655E-01, -2.12998316E-01, 6.84030479E-02,
        6.84030479E-02, 7.50192517E-01,
    ])
    assert abs(mol.wfn.exp_alpha.energies - expected_energies).max() < 3e-5

    ham.compute()
    # compare with g09
    assert abs(ham.cache['energy_ne'] - -2.981579553570E+02) < 1e-5
    assert abs(ham.cache['energy_kin'] - 1.061620887711E+02) < 1e-5
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_x_dirac'] - 6.247259253877E+01) < 1e-4
    assert abs(ham.cache['energy'] - -106.205213597) < 1e-4
    assert abs(ham.cache['energy_nn'] - 23.3180604505) < 1e-8


def test_fock_h3_hfs_321g():
    # The fock operator is tested by doing an SCF and checking the converged
    # energies
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    mol = Molecule.from_file(fn_fchk)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'veryfine', random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        UnrestrictedOneBodyTerm(kin, 'kin'),
        UnrestrictedDirectTerm(er, 'hartree'),
        UnrestrictedGridGroup(mol.obasis, grid, [
            UnrestrictedDiracExchange(),
        ]),
        UnrestrictedOneBodyTerm(na, 'ne'),
    ]
    external = {'nn': compute_nucnuc(mol.coordinates, mol.pseudo_numbers)}
    ham = UnrestrictedEffectiveHamiltonian(terms, external)

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 1e-6

    # Converge from scratch
    guess_core_hamiltonian(mol.wfn, olp, kin, na)
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) > 1e-8
    converge_scf(ham, mol.wfn, mol.lf, olp)
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 1e-8

    # test orbital energies
    expected_energies = np.array([
        -4.93959157E-01, -1.13961330E-01, 2.38730924E-01, 7.44216538E-01,
        8.30143356E-01, 1.46613581E+00
    ])
    assert abs(mol.wfn.exp_alpha.energies - expected_energies).max() < 1e-5
    expected_energies = np.array([
        -4.34824166E-01, 1.84114514E-04, 3.24300545E-01, 7.87622756E-01,
        9.42415831E-01, 1.55175481E+00
    ])
    assert abs(mol.wfn.exp_beta.energies - expected_energies).max() < 1e-5

    ham.compute()
    # compare with g09
    assert abs(ham.cache['energy_ne'] - -6.832069993374E+00) < 1e-5
    assert abs(ham.cache['energy_kin'] - 1.870784279014E+00) < 1e-5
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_x_dirac'] - 1.658810998195E+00) < 1e-6
    assert abs(ham.cache['energy'] - -1.412556114057104E+00) < 1e-5
    assert abs(ham.cache['energy_nn'] - 1.8899186021) < 1e-8


def test_cubic_interpolation_hfs_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    mol = Molecule.from_file(fn_fchk)

    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, random_rotate=False)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RestrictedOneBodyTerm(kin, 'kin'),
        RestrictedDirectTerm(er, 'hartree'),
        RestrictedGridGroup(mol.obasis, grid, [
            RestrictedDiracExchange(),
        ]),
        RestrictedOneBodyTerm(na, 'ne'),
    ]
    ham = RestrictedEffectiveHamiltonian(terms)

    dm0 = mol.lf.create_one_body()
    dm0.assign(mol.wfn.dm_alpha)
    guess_core_hamiltonian(mol.wfn, olp, kin, na)
    dm1 = mol.lf.create_one_body()
    dm1.assign(mol.wfn.dm_alpha)

    check_cubic_cs_wrapper(ham, dm0, dm1)


def test_perturbation():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)

    # Without perturbation
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RestrictedOneBodyTerm(kin, 'kin'),
        RestrictedDirectTerm(er, 'hartree'),
        RestrictedExchangeTerm(er, 'x_hf'),
        RestrictedOneBodyTerm(na, 'ne'),
    ]
    ham = RestrictedEffectiveHamiltonian(terms)
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) > 1e-8
    converge_scf(ham, mol.wfn, mol.lf, olp)
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 1e-8
    energy0 = ham.compute()

    # Construct a perturbation based on the Mulliken AIM operator
    assert mol.obasis.nbasis % 2 == 0
    nfirst = mol.obasis.nbasis / 2
    operator = mol.obasis.compute_overlap(mol.lf).copy()
    operator._array[:nfirst,nfirst:] *= 0.5
    operator._array[nfirst:,:nfirst] *= 0.5
    operator._array[nfirst:,nfirst:] = 0.0

    # Apply the perturbation with oposite signs and check that, because of
    # symmetry, the energy of the perturbed wavefunction is the same in both
    # cases, and higher than the unperturbed.
    energy1_old = None
    for scale in 0.1, -0.1:
        # Perturbation
        tmp = operator.copy()
        tmp.iscale(scale)
        perturbation = RestrictedOneBodyTerm(tmp, 'pert')
        # Hamiltonian
        terms = [
            RestrictedOneBodyTerm(kin, 'kin'),
            RestrictedDirectTerm(er, 'hartree'),
            RestrictedExchangeTerm(er, 'x_hf'),
            RestrictedOneBodyTerm(na, 'ne'),
            perturbation,
        ]
        ham = RestrictedEffectiveHamiltonian(terms)
        assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) > 1e-8
        converge_scf_oda(ham, mol.wfn, mol.lf, olp)
        assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 1e-8
        energy1 = ham.compute()
        energy1 -= ham.cache['energy_pert']

        assert energy1 > energy0
        if energy1_old is None:
            energy1_old = energy1
        else:
            assert abs(energy1 - energy1_old) < 1e-7


def test_ghost_hf():
    fn_fchk = context.get_fn('test/water_dimer_ghost.fchk')
    mol = Molecule.from_file(fn_fchk)
    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        RestrictedOneBodyTerm(kin, 'kin'),
        RestrictedDirectTerm(er, 'hartree'),
        RestrictedExchangeTerm(er, 'x_hf'),
        RestrictedOneBodyTerm(na, 'ne'),
    ]
    ham = RestrictedEffectiveHamiltonian(terms)
    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 1e-5
