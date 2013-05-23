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


import numpy as np
from horton import *
from horton.meanfield.test.common import *


def test_hamiltonian_init():
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    numbers = np.array([1, 1])
    sys = System(coordinates, numbers, obasis='STO-3G')
    sys.init_wfn(0, 1)

    # test if no terms gives a ValueError
    try:
        ham = Hamiltonian(sys, [])
        assert False
    except ValueError:
        pass

    # test if terms are added automagically
    ham = Hamiltonian(sys, [HartreeFock()])
    assert sum(isinstance(term, KineticEnergy) for term in ham.terms) == 1
    assert sum(isinstance(term, ExternalPotential) for term in ham.terms) == 1
    assert sum(isinstance(term, Hartree) for term in ham.terms) == 1

    # test if the necessary operators are constructed in the system object
    assert 'kin' in sys.operators
    assert 'na' in sys.operators
    assert 'er' in sys.operators
    assert 'olp' in sys.operators

    # check attribute of HartreeFock term
    assert ham.terms[0].fraction_exchange == 1.0


def test_energy_hydrogen():
    fn_fchk = context.get_fn('test/h_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    ham = Hamiltonian(sys, [HartreeFock()])
    ham.compute_energy()
    assert abs(sys.props['energy'] - -4.665818503844346E-01) < 1e-8


def test_energy_n2_hfs_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), DiracExchange()], grid)
    ham.compute_energy()

    # Compare energies
    assert abs(sys.props['energy_ne'] - -2.981579553570E+02) < 1e-6
    assert abs(sys.props['energy_kin'] - 1.061620887711E+02) < 1e-6
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_dirac'] - 6.247259253877E+01) < 1e-4
    assert abs(sys.props['energy'] - -106.205213597) < 1e-4
    assert abs(sys.props['energy_nn'] - 23.3180604505) < 3e-9


    # Test if the grid potential data is properly converted into an operator:
    ev1 = grid.integrate(ham.cache.load('pot_exchange_dirac_alpha'), ham.cache.load('rho_alpha'))
    dma = sys.lf.create_one_body()
    ev2 = ham.cache.load('op_exchange_dirac_alpha').expectation_value(sys.wfn.dm_alpha)
    assert abs(ev1 - ev2) < 1e-10

    # When repeating, we should get the same
    ham.compute_energy()
    assert abs(sys.props['energy'] - -106.205213597) < 1e-4

    # check symmetry
    ham.cache.load('op_exchange_dirac_alpha').check_symmetry()


def test_fock_n2_hfs_sto3g():
    # The fock operator is tested by doing an SCF an checking the converged
    # energies
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), DiracExchange()], grid)

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error(ham) < 1e-5

    # Converge from scratch
    guess_hamiltonian_core(sys)
    assert convergence_error(ham) > 1e-8
    assert converge_scf(ham)
    assert convergence_error(ham) < 1e-8

    # test orbital energies
    expected_energies = np.array([
        -1.37107053E+01, -1.37098006E+01, -9.60673085E-01, -3.57928483E-01,
        -3.16017655E-01, -3.16017655E-01, -2.12998316E-01, 6.84030479E-02,
        6.84030479E-02, 7.50192517E-01,
    ])
    assert abs(sys.wfn.exp_alpha.energies - expected_energies).max() < 3e-5

    ham.compute_energy()
    # compare with g09
    assert abs(sys.props['energy_ne'] - -2.981579553570E+02) < 1e-5
    assert abs(sys.props['energy_kin'] - 1.061620887711E+02) < 1e-5
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_dirac'] - 6.247259253877E+01) < 1e-4
    assert abs(sys.props['energy'] - -106.205213597) < 1e-4
    assert abs(sys.props['energy_nn'] - 23.3180604505) < 1e-8


def test_fock_h3_hfs_321g():
    # The fock operator is tested by doing an SCF an checking the converged
    # energies
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)
    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), DiracExchange()], grid)

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error(ham) < 1e-6

    # Converge from scratch
    guess_hamiltonian_core(sys)
    assert convergence_error(ham) > 1e-8
    assert converge_scf(ham)
    assert convergence_error(ham) < 1e-8

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

    ham.compute_energy()
    # compare with g09
    assert abs(sys.props['energy_ne'] - -6.832069993374E+00) < 1e-5
    assert abs(sys.props['energy_kin'] - 1.870784279014E+00) < 1e-5
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_dirac'] - 1.658810998195E+00) < 1e-6
    assert abs(sys.props['energy'] - -1.412556114057104E+00) < 1e-5
    assert abs(sys.props['energy_nn'] - 1.8899186021) < 1e-8


def test_cubic_interpolation_hfs_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), DiracExchange()], grid)

    dm0 = sys.lf.create_one_body()
    dm0.assign(sys.wfn.dm_alpha)
    guess_hamiltonian_core(sys)
    dm1 = sys.lf.create_one_body()
    dm1.assign(sys.wfn.dm_alpha)

    check_cubic_cs_wrapper(ham, dm0, dm1)


def test_external_potential_copy():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)
    ep = ExternalPotential()
    ep.prepare_system(sys, None, None)
    assert not (ep.operator is sys.get_nuclear_attraction())


def test_custom_term():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)

    # Without perturbation
    ham = Hamiltonian(sys, [HartreeFock()])
    assert convergence_error(ham) > 1e-8
    assert converge_scf(ham)
    assert convergence_error(ham) < 1e-8
    energy0 = ham.compute_energy()

    # Construct a perturbation baed on the Mulliken AIM operator
    assert sys.obasis.nbasis % 2 == 0
    nfirst = sys.obasis.nbasis / 2
    operator = ham.overlap.copy()
    operator._array[:nfirst,nfirst:] *= 0.5
    operator._array[nfirst:,:nfirst] *= 0.5
    operator._array[nfirst:,nfirst:] = 0.0

    # Apply the perturbation with oposite signs and check that, because of
    # symmetry, the energy of the perturbed wavefunction is the same in both
    # cases, and higher than the unperturbed.
    energy1_old = None
    for scale in 0.1, -0.1:
        # With perturbation
        tmp = operator.copy()
        tmp.iscale(scale)
        perturbation = CustomFixedTerm(tmp, 'pert')
        ham = Hamiltonian(sys, [HartreeFock(), perturbation])
        assert convergence_error(ham) > 1e-8
        assert converge_scf_oda(ham)
        assert convergence_error(ham) < 1e-8
        energy1 = ham.compute_energy()
        energy1 -= sys.props['energy_pert']

        assert energy1 > energy0
        if energy1_old is None:
            energy1_old = energy1
        else:
            assert abs(energy1 - energy1_old) < 1e-7


def test_auto_complete():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)

    # HF case
    ham = Hamiltonian(sys, [HartreeFock()])
    assert any(isinstance(term, KineticEnergy) for term in ham.terms)
    assert any(isinstance(term, ExternalPotential) for term in ham.terms)

    # DFT case
    grid = BeckeMolGrid(sys)
    ham = Hamiltonian(sys, [DiracExchange()], grid)
    assert any(isinstance(term, KineticEnergy) for term in ham.terms)
    assert any(isinstance(term, ExternalPotential) for term in ham.terms)
    assert any(isinstance(term, Hartree) for term in ham.terms)

    # special behavior
    ham = Hamiltonian(sys, [HartreeFock()], auto_complete=False)
    assert not any(isinstance(term, KineticEnergy) for term in ham.terms)
    assert not any(isinstance(term, ExternalPotential) for term in ham.terms)
