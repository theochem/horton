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
    sys = System.from_file(fn_fchk)
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
        DirectTerm(er, sys.lf, sys.wfn),
        ExchangeTerm(er, sys.lf, sys.wfn),
        OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
    ]
    external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
    ham = Hamiltonian(sys, terms, external)
    ham.compute()
    assert abs(ham.cache['energy'] - -4.665818503844346E-01) < 1e-8


def test_energy_n2_hfs_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers, random_rotate=False)
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
        DirectTerm(er, sys.lf, sys.wfn),
        GridGroup(sys.obasis, grid, sys.lf, sys.wfn, [
            DiracExchange(sys.wfn),
        ]),
        OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
    ]
    external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
    ham = Hamiltonian(sys, terms, external)
    ham.compute()

    # Compare energies
    assert abs(ham.cache['energy_ne'] - -2.981579553570E+02) < 1e-6
    assert abs(ham.cache['energy_kin'] - 1.061620887711E+02) < 1e-6
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_exchange_dirac'] - 6.247259253877E+01) < 1e-4
    assert abs(ham.cache['energy'] - -106.205213597) < 1e-4
    assert abs(ham.cache['energy_nn'] - 23.3180604505) < 3e-9


    # Test if the grid potential data can be properly converted into an operator:
    pot = ham.cache.load('pot_exchange_dirac_alpha')
    ev1 = grid.integrate(pot, ham.cache.load('rho_alpha'))
    op = sys.lf.create_one_body()
    sys.compute_grid_density_fock(grid.points, grid.weights, pot, op)
    ev2 = op.expectation_value(sys.wfn.dm_alpha)
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
    sys = System.from_file(fn_fchk)
    grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers, 'veryfine', random_rotate=False)
    olp = sys.get_overlap()
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
        DirectTerm(er, sys.lf, sys.wfn),
        GridGroup(sys.obasis, grid, sys.lf, sys.wfn, [
            DiracExchange(sys.wfn),
        ]),
        OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
    ]
    external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
    ham = Hamiltonian(sys, terms, external)

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham) < 1e-5

    # Converge from scratch
    guess_core_hamiltonian(sys.wfn, olp, kin, nai)
    assert convergence_error_eigen(ham) > 1e-8
    converge_scf(ham)
    assert convergence_error_eigen(ham) < 1e-8

    # test orbital energies
    expected_energies = np.array([
        -1.37107053E+01, -1.37098006E+01, -9.60673085E-01, -3.57928483E-01,
        -3.16017655E-01, -3.16017655E-01, -2.12998316E-01, 6.84030479E-02,
        6.84030479E-02, 7.50192517E-01,
    ])
    assert abs(sys.wfn.exp_alpha.energies - expected_energies).max() < 3e-5

    ham.compute()
    # compare with g09
    assert abs(ham.cache['energy_ne'] - -2.981579553570E+02) < 1e-5
    assert abs(ham.cache['energy_kin'] - 1.061620887711E+02) < 1e-5
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_exchange_dirac'] - 6.247259253877E+01) < 1e-4
    assert abs(ham.cache['energy'] - -106.205213597) < 1e-4
    assert abs(ham.cache['energy_nn'] - 23.3180604505) < 1e-8


def test_fock_h3_hfs_321g():
    # The fock operator is tested by doing an SCF and checking the converged
    # energies
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)
    grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers, 'veryfine', random_rotate=False)
    olp = sys.get_overlap()
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
        DirectTerm(er, sys.lf, sys.wfn),
        GridGroup(sys.obasis, grid, sys.lf, sys.wfn, [
            DiracExchange(sys.wfn),
        ]),
        OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
    ]
    external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
    ham = Hamiltonian(sys, terms, external)

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham) < 1e-6

    # Converge from scratch
    guess_core_hamiltonian(sys.wfn, olp, kin, nai)
    assert convergence_error_eigen(ham) > 1e-8
    converge_scf(ham)
    assert convergence_error_eigen(ham) < 1e-8

    # test orbital energies
    expected_energies = np.array([
        -4.93959157E-01, -1.13961330E-01, 2.38730924E-01, 7.44216538E-01,
        8.30143356E-01, 1.46613581E+00
    ])
    assert abs(sys.wfn.exp_alpha.energies - expected_energies).max() < 1e-5
    expected_energies = np.array([
        -4.34824166E-01, 1.84114514E-04, 3.24300545E-01, 7.87622756E-01,
        9.42415831E-01, 1.55175481E+00
    ])
    assert abs(sys.wfn.exp_beta.energies - expected_energies).max() < 1e-5

    ham.compute()
    # compare with g09
    assert abs(ham.cache['energy_ne'] - -6.832069993374E+00) < 1e-5
    assert abs(ham.cache['energy_kin'] - 1.870784279014E+00) < 1e-5
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_exchange_dirac'] - 1.658810998195E+00) < 1e-6
    assert abs(ham.cache['energy'] - -1.412556114057104E+00) < 1e-5
    assert abs(ham.cache['energy_nn'] - 1.8899186021) < 1e-8


def test_cubic_interpolation_hfs_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)

    grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers, random_rotate=False)
    olp = sys.get_overlap()
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
        DirectTerm(er, sys.lf, sys.wfn),
        GridGroup(sys.obasis, grid, sys.lf, sys.wfn, [
            DiracExchange(sys.wfn),
        ]),
        OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
    ]
    ham = Hamiltonian(sys, terms)

    dm0 = sys.lf.create_one_body()
    dm0.assign(sys.wfn.dm_alpha)
    guess_core_hamiltonian(sys.wfn, olp, kin, nai)
    dm1 = sys.lf.create_one_body()
    dm1.assign(sys.wfn.dm_alpha)

    check_cubic_cs_wrapper(ham, dm0, dm1)


def test_perturbation():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)

    # Without perturbation
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
        DirectTerm(er, sys.lf, sys.wfn),
        ExchangeTerm(er, sys.lf, sys.wfn),
        OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
    ]
    ham = Hamiltonian(sys, terms)
    assert convergence_error_eigen(ham) > 1e-8
    converge_scf(ham)
    assert convergence_error_eigen(ham) < 1e-8
    energy0 = ham.compute()

    # Construct a perturbation based on the Mulliken AIM operator
    assert sys.obasis.nbasis % 2 == 0
    nfirst = sys.obasis.nbasis / 2
    operator = sys.get_overlap().copy()
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
        perturbation = OneBodyTerm(tmp, sys.lf, sys.wfn, 'pert')
        # Hamiltonian
        terms = [
            OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
            DirectTerm(er, sys.lf, sys.wfn),
            ExchangeTerm(er, sys.lf, sys.wfn),
            OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
            perturbation,
        ]
        ham = Hamiltonian(sys, terms)
        assert convergence_error_eigen(ham) > 1e-8
        converge_scf_oda(ham)
        assert convergence_error_eigen(ham) < 1e-8
        energy1 = ham.compute()
        energy1 -= ham.cache['energy_pert']

        assert energy1 > energy0
        if energy1_old is None:
            energy1_old = energy1
        else:
            assert abs(energy1 - energy1_old) < 1e-7


def test_add_term():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)
    ham = Hamiltonian(sys, [ExchangeTerm(sys.get_electron_repulsion(),
                                         sys.lf, sys.wfn)])
    term = OneBodyTerm(sys.get_kinetic(), sys.lf, sys.wfn, 'kin')
    ham.add_term(term)
    assert term._hamiltonian is ham


def test_ghost_hf():
    fn_fchk = context.get_fn('test/water_dimer_ghost.fchk')
    sys = System.from_file(fn_fchk)
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, sys.lf, sys.wfn, 'kin'),
        DirectTerm(er, sys.lf, sys.wfn),
        ExchangeTerm(er, sys.lf, sys.wfn),
        OneBodyTerm(nai, sys.lf, sys.wfn, 'ne'),
    ]
    ham = Hamiltonian(sys, terms)
    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error_eigen(ham) < 1e-5
