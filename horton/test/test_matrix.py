# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


def get_water_sto3g_hf(lf=None):
    if lf is None:
        lf = DenseLinalgFactory()
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
    overlap, kinetic, nuclear_attraction, electronic_repulsion = load_operators_g09(fn, lf)
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
    wfn = ClosedShellWFN(nep=5, lf=lf, nbasis=7)
    wfn.expansion.coeffs[:] = coeffs
    wfn.expansion.energies[:] = epsilons
    return lf, overlap, kinetic, nuclear_attraction, electronic_repulsion, wfn


def test_fock_matrix_eigen():
    lf, overlap, kinetic, nuclear_attraction, electronic_repulsion, wfn = get_water_sto3g_hf()
    nbasis = overlap.nbasis

    dm = lf.create_one_body(7)
    coulomb = lf.create_one_body(nbasis)
    exchange = lf.create_one_body(nbasis)
    wfn.compute_alpha_density_matrix(dm)
    electronic_repulsion.apply_direct(dm, coulomb)
    electronic_repulsion.apply_exchange(dm, exchange)

    # Construct the Fock operator
    fock = lf.create_one_body(nbasis)
    fock.iadd(kinetic, 1)
    fock.iadd(nuclear_attraction, -1)
    fock.iadd(coulomb, 2)
    fock.iadd(exchange, -1)

    # Check for convergence
    expansion = wfn.expansion
    error = lf.error_eigen(fock, overlap, expansion)
    assert error > 0
    assert error < 1e-4

    # Check self-consistency of the orbital energies
    old_energies = expansion.energies.copy()
    lf.diagonalize(fock, overlap, expansion)
    assert abs(expansion.energies - old_energies).max() < 1e-4



def test_kinetic_energy_water_sto3g():
    lf, overlap, kinetic, nuclear_attraction, electronic_repulsion, wfn = get_water_sto3g_hf()
    dm = lf.create_one_body(7)
    wfn.compute_density_matrix(dm)
    ekin = kinetic.expectation_value(dm)
    assert abs(ekin - 74.60736832935) < 1e-4


def test_ortho_water_sto3g():
    lf, overlap, kinetic, nuclear_attraction, electronic_repulsion, wfn = get_water_sto3g_hf()
    for i0 in xrange(7):
        orb0 = wfn.expansion.coeffs[:,i0]
        for i1 in xrange(i0+1):
            orb1 = wfn.expansion.coeffs[:,i1]
            check = overlap.dot(orb0, orb1)
            assert abs(check - (i0==i1)) < 1e-4


def test_potential_energy_water_sto3g_hf():
    lf, overlap, kinetic, nuclear_attraction, electronic_repulsion, wfn = get_water_sto3g_hf()
    dm = lf.create_one_body(7)
    wfn.compute_density_matrix(dm)
    epot = -nuclear_attraction.expectation_value(dm)
    assert abs(epot - (-197.1170963957)) < 2e-3


def test_electron_electron_water_sto3g_hf():
    lf, overlap, kinetic, nuclear_attraction, electronic_repulsion, wfn = get_water_sto3g_hf()
    dm = lf.create_one_body(7)
    coulomb = lf.create_one_body(7)
    exchange = lf.create_one_body(7)
    wfn.compute_alpha_density_matrix(dm)
    electronic_repulsion.apply_direct(dm, coulomb)
    electronic_repulsion.apply_exchange(dm, exchange)
    eee = 2*coulomb.expectation_value(dm) \
          - exchange.expectation_value(dm)
    assert abs(eee - 38.29686853319) < 1e-4


def test_hartree_fock_water():
    lf, overlap, kinetic, nuclear_attraction, electronic_repulsion, wfn0 = get_water_sto3g_hf()
    nbasis = overlap.nbasis

    # Construct a wavefunction
    wfn = ClosedShellWFN(nep=5, lf=lf, nbasis=nbasis)

    # Construct the hamiltonian core guess
    hamcore = lf.create_one_body(nbasis)
    hamcore.iadd(kinetic, 1)
    hamcore.iadd(nuclear_attraction, -1)
    lf.diagonalize(hamcore, overlap, wfn.expansion)

    # The SCF loop
    coulomb = lf.create_one_body(nbasis)
    exchange = lf.create_one_body(nbasis)
    fock = lf.create_one_body(nbasis)
    dm = lf.create_one_body(nbasis)
    for i in xrange(1000):
        # Construct the Fock operator
        fock.reset()
        fock.iadd(hamcore, 1)
        wfn.compute_alpha_density_matrix(dm)
        electronic_repulsion.apply_direct(dm, coulomb)
        electronic_repulsion.apply_exchange(dm, exchange)
        fock.iadd(coulomb, 2)
        fock.iadd(exchange, -1)
        # Check for convergence
        error = lf.error_eigen(fock, overlap, wfn.expansion)
        if error < 1e-10:
            break
        # Diagonalize the fock operator
        lf.diagonalize(fock, overlap, wfn.expansion)

    assert abs(wfn.expansion.energies - wfn0.expansion.energies).max() < 1e-4

    # Check the hartree-fock energy
    hf1 = sum([
        -2*coulomb.expectation_value(dm),
        +1*exchange.expectation_value(dm),
    ]) + wfn.expansion.energies[:wfn.nep].sum()*2
    hf2 = sum([
        2*kinetic.expectation_value(dm),
        -2*nuclear_attraction.expectation_value(dm),
        +2*coulomb.expectation_value(dm),
        -exchange.expectation_value(dm),
    ])
    enn = 9.2535672047 # nucleus-nucleus interaction
    assert abs(hf1 + enn - (-74.9592923284)) < 1e-4
    assert abs(hf2 + enn - (-74.9592923284)) < 1e-4


def test_dense_one_body_trace():
    lf = DenseLinalgFactory()
    op1 = lf.create_one_body(3)
    op1._array[:] = np.random.uniform(-1, 1, (3,3))
    op1._valid = True
    assert op1.trace() == op1._array[0,0] + op1._array[1,1] + op1._array[2,2]


def test_dense_one_body_itranspose():
    lf = DenseLinalgFactory()
    op1 = lf.create_one_body(3)
    op2 = lf.create_one_body(3)
    op1._array[:] = np.random.uniform(-1, 1, (3,3))
    op2._array[:] = op1._array
    op2.itranspose()
    assert op1._array[0,1] == op2._array[1,0]


def test_dense_one_body_iscale():
    lf = DenseLinalgFactory()
    op = lf.create_one_body(3)
    op._array[:] = np.random.uniform(-1, 1, (3,3))
    tmp = op._array.copy()
    op.iscale(3.0)
    assert abs(op._array - 3*tmp).max() < 1e-10


def test_dense_expansion_properties():
    lf = DenseLinalgFactory()
    ex = lf.create_expansion(10, 8, do_energies=False)
    assert ex.nbasis == 10
    assert ex.nfn == 8
    assert ex.coeffs.shape == (10,8) # orbitals stored as columns
    assert ex.energies is None


def test_dense_one_body_properties():
    lf = DenseLinalgFactory()
    op = lf.create_one_body(3)
    assert op.nbasis == 3
    assert not op.valid
    op.set_element(0, 1, 1.2)
    assert op.valid
    op.invalidate()
    assert not op.valid


def test_dense_two_body_properties():
    lf = DenseLinalgFactory()
    op = lf.create_two_body(3)
    assert op.nbasis == 3
