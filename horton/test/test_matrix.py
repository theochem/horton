# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011 Toon Verstraelen <Toon.Verstraelen@UGent.be>, ...
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


from horton import *


def test_hartree_fock_water():
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    overlap, kinetic, nuclear_attraction, electron_repulsion = load_operators_g09(fn)
    nbasis = overlap.size

    # Construct a wavefunction
    wfn = ClosedShellWFN(nep=5, nbasis=nbasis)

    # Construct the hamiltonian core guess
    hamcore = DenseOneBody(nbasis)
    hamcore.iadd(kinetic, 1)
    hamcore.iadd(nuclear_attraction, -1)
    diagonalize(hamcore, overlap, wfn)

    # The SCF loop
    coulomb = DenseOneBody(nbasis)
    exchange = DenseOneBody(nbasis)
    fock = DenseOneBody(nbasis)
    for i in xrange(1000):
        # Construct the Fock operator
        fock.reset()
        fock.iadd(hamcore, 1)
        wfn.apply_two_body(electron_repulsion, coulomb, exchange)
        fock.iadd(coulomb, 2)
        fock.iadd(exchange, -1)
        # Check for convergence
        error = error_eigen(fock, overlap, wfn)
        if error < 1e-10:
            break
        # Diagonalize the fock operator
        diagonalize(fock, overlap, wfn)

    check_epsilons = np.array([
        -2.02333942E+01, -1.26583942E+00, -6.29365088E-01, -4.41724988E-01,
        -3.87671783E-01, 6.03082408E-01, 7.66134805E-01
    ])
    assert abs(wfn._epsilons - check_epsilons).max() < 1e-4

    # Check the hartree-fock energy
    dm = wfn.get_density_matrix()
    wfn.apply_two_body(electron_repulsion, coulomb, exchange)
    hf1 = sum([
        -2*coulomb.expectation_value(dm),
        +1*exchange.expectation_value(dm),
    ]) + wfn._epsilons[:wfn.nep].sum()*2
    hf2 = sum([
        2*kinetic.expectation_value(dm),
        -2*nuclear_attraction.expectation_value(dm),
        +2*coulomb.expectation_value(dm),
        -exchange.expectation_value(dm),
    ])
    enn = 9.2535672047 # nucleus-nucleus interaction
    assert abs(hf1 + enn - (-74.9592923284)) < 1e-4
    assert abs(hf2 + enn - (-74.9592923284)) < 1e-4


def get_water_sto3g_hf_wfn():
    coeffs = np.array([
        9.94099882E-01, 2.67799213E-02, 3.46630004E-03, -1.54676269E-15,
        2.45105601E-03, -6.08393842E-03, -6.08393693E-03, -2.32889095E-01,
        8.31788042E-01, 1.03349385E-01, 9.97532839E-17, 7.30794097E-02,
        1.60223990E-01, 1.60223948E-01, 1.65502862E-08, -9.03020258E-08,
        -3.46565859E-01, -2.28559667E-16, 4.90116062E-01, 4.41542336E-01,
        -4.41542341E-01, 1.00235366E-01, -5.23423149E-01, 6.48259144E-01,
        -5.78009326E-16, 4.58390414E-01, 2.69085788E-01, 2.69085849E-01,
        8.92936017E-17, -1.75482465E-16, 2.47517845E-16, 1.00000000E+00,
        5.97439610E-16, -3.70474007E-17, -2.27323914E-17, -1.35631600E-01,
        9.08581133E-01, 5.83295647E-01, -4.37819173E-16, 4.12453695E-01,
        -8.07337352E-01, -8.07337875E-01, 5.67656309E-08, -4.29452066E-07,
        5.82525068E-01, -6.76605679E-17, -8.23811720E-01, 8.42614916E-01,
        -8.42614243E-01
    ]).reshape(7,7).T
    epsilons = np.array([
        -2.02333942E+01, -1.26583942E+00, -6.29365088E-01, -4.41724988E-01,
        -3.87671783E-01, 6.03082408E-01, 7.66134805E-01
    ])
    wfn = ClosedShellWFN(nep=5, nbasis=7)
    wfn._expansion._coeffs[:] = coeffs
    wfn._epsilons[:] = epsilons
    return wfn


def test_fock_matrix_eigen():
    wfn = get_water_sto3g_hf_wfn()
    dm = wfn.get_density_matrix()
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    overlap, kinetic, nuclear_attraction, electron_repulsion = load_operators_g09(fn)
    nbasis = overlap.size

    coulomb = DenseOneBody(nbasis)
    exchange = DenseOneBody(nbasis)
    wfn.apply_two_body(electron_repulsion, coulomb, exchange)

    # Construct the Fock operator
    fock = DenseOneBody(nbasis)
    fock.iadd(kinetic, 1)
    fock.iadd(nuclear_attraction, -1)
    fock.iadd(coulomb, 2)
    fock.iadd(exchange, -1)

    # Check for convergence
    error = error_eigen(fock, overlap, wfn)
    assert error > 0
    assert error < 1e-4

    # Check self-consistency of the orbital energies
    old_epsilons = wfn._epsilons.copy()
    diagonalize(fock, overlap, wfn)
    assert abs(wfn._epsilons - old_epsilons).max() < 1e-4



def test_kinetic_energy_water_sto3g():
    wfn = get_water_sto3g_hf_wfn()
    dm = wfn.get_density_matrix()
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    overlap, kinetic, nuclear_attraction, electron_repulsion = load_operators_g09(fn)
    ekin = 2*kinetic.expectation_value(dm)
    assert abs(ekin - 74.60736832935) < 1e-4


def test_ortho_water_sto3g():
    wfn = get_water_sto3g_hf_wfn()
    dm = wfn.get_density_matrix()
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    overlap, kinetic, nuclear_attraction, electron_repulsion = load_operators_g09(fn)
    for i0 in xrange(7):
        orb0 = wfn._expansion._coeffs[:,i0]
        for i1 in xrange(i0+1):
            orb1 = wfn._expansion._coeffs[:,i1]
            check = np.dot(orb0, np.dot(overlap._array, orb1))
            assert abs(check - (i0==i1)) < 1e-4


def test_potential_energy_water_sto3g_hf():
    wfn = get_water_sto3g_hf_wfn()
    dm = wfn.get_density_matrix()
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    overlap, kinetic, nuclear_attraction, electron_repulsion = load_operators_g09(fn)
    epot = -2*nuclear_attraction.expectation_value(dm)
    assert abs(epot - (-197.1170963957)) < 2e-3


def test_electron_electron_water_sto3g_hf():
    wfn = get_water_sto3g_hf_wfn()
    dm = wfn.get_density_matrix()
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    overlap, kinetic, nuclear_attraction, electron_repulsion = load_operators_g09(fn)
    coulomb = DenseOneBody(7)
    exchange = DenseOneBody(7)
    wfn.apply_two_body(electron_repulsion, coulomb, exchange)
    eee = 2*coulomb.expectation_value(dm) - exchange.expectation_value(dm)
    assert abs(eee - 38.29686853319) < 1e-4
