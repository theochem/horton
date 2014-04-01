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

from horton import *


def get_water_sto3g_hf(lf=None):
    if lf is None:
        lf = DenseLinalgFactory()
    fn = context.get_fn('test/water_sto3g_hf_g03.log')
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
    occ_model = AufbauOccModel(5)
    wfn = RestrictedWFN(lf, 7, occ_model)
    wfn.init_exp('alpha')
    exp_alpha = wfn.exp_alpha
    exp_alpha.coeffs[:] = coeffs
    exp_alpha.energies[:] = epsilons
    occ_model.assign(exp_alpha)
    assert (exp_alpha.occupations == np.array([1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])).all()
    # convert the cache dictionary to a real cache object
    cache = Cache()
    data = load_operators_g09(fn, lf)
    for key, value in data.iteritems():
        cache.dump(key, value)
    return lf, cache, wfn


def test_fock_matrix_eigen():
    lf, cache, wfn = get_water_sto3g_hf()
    nbasis = cache['olp'].nbasis

    hartree = lf.create_one_body(nbasis)
    exchange = lf.create_one_body(nbasis)
    dm = wfn.dm_alpha
    cache['er'].apply_direct(dm, hartree)
    cache['er'].apply_exchange(dm, exchange)

    # Construct the Fock operator
    fock = lf.create_one_body(nbasis)
    fock.iadd(cache['kin'], 1)
    fock.iadd(cache['na'], -1)
    fock.iadd(hartree, 2)
    fock.iadd(exchange, -1)

    # Check for convergence
    exp_alpha = wfn.exp_alpha
    error = lf.error_eigen(fock, cache['olp'], exp_alpha)
    assert error > 0
    assert error < 1e-4

    # Check self-consistency of the orbital energies
    old_energies = exp_alpha.energies.copy()
    exp_alpha.derive_from_fock_matrix(fock, cache['olp'])
    assert abs(exp_alpha.energies - old_energies).max() < 1e-4



def test_kinetic_energy_water_sto3g():
    lf, cache, wfn = get_water_sto3g_hf()
    dm = wfn.dm_full
    ekin = cache['kin'].expectation_value(dm)
    assert abs(ekin - 74.60736832935) < 1e-4


def test_ortho_water_sto3g():
    lf, cache, wfn = get_water_sto3g_hf()
    exp_alpha = wfn.exp_alpha
    for i0 in xrange(7):
        orb0 = exp_alpha.coeffs[:,i0]
        for i1 in xrange(i0+1):
            orb1 = exp_alpha.coeffs[:,i1]
            check = cache['olp'].dot(orb0, orb1)
            assert abs(check - (i0==i1)) < 1e-4


def test_potential_energy_water_sto3g_hf():
    lf, cache, wfn = get_water_sto3g_hf()
    dm = wfn.dm_full
    #epot = -nuclear_attraction.expectation_value(dm)
    epot = -cache['na'].expectation_value(dm)
    assert abs(epot - (-197.1170963957)) < 2e-3


def test_electron_electron_water_sto3g_hf():
    lf, cache, wfn = get_water_sto3g_hf()
    hartree = lf.create_one_body(7)
    exchange = lf.create_one_body(7)
    dm = wfn.dm_alpha
    cache['er'].apply_direct(dm, hartree)
    cache['er'].apply_exchange(dm, exchange)
    eee = 2*hartree.expectation_value(dm) \
          - exchange.expectation_value(dm)
    assert abs(eee - 38.29686853319) < 1e-4


def test_hartree_fock_water():
    lf, cache, wfn0 = get_water_sto3g_hf()
    nbasis = cache['olp'].nbasis

    # Construct a wavefunction
    occ_model = AufbauOccModel(5)
    wfn = RestrictedWFN(lf, nbasis, occ_model)

    # Construct the hamiltonian core guess
    hamcore = lf.create_one_body(nbasis)
    hamcore.iadd(cache['kin'], 1)
    hamcore.iadd(cache['na'], -1)
    wfn.clear()
    exp_alpha1 = wfn.update_exp(hamcore, cache['olp'])
    assert (exp_alpha1.energies != 0.0).any()


    # The SCF loop
    hartree = lf.create_one_body(nbasis)
    exchange = lf.create_one_body(nbasis)
    fock = lf.create_one_body(nbasis)
    #dm = lf.create_one_body(nbasis)
    for i in xrange(1000):
        # Construct the Fock operator
        fock.clear()
        fock.iadd(hamcore, 1)
        cache['er'].apply_direct(wfn.dm_alpha, hartree)
        cache['er'].apply_exchange(wfn.dm_alpha, exchange)
        fock.iadd(hartree, 2)
        fock.iadd(exchange, -1)
        # Check for convergence
        error = lf.error_eigen(fock, cache['olp'], wfn.exp_alpha)
        if error < 1e-10:
            break
        # Derive the expansion and the density matrix from the fock operator
        wfn.clear()
        wfn.update_exp(fock, cache['olp'])

    exp_alpha = wfn.exp_alpha
    exp_alpha0 = wfn0.exp_alpha
    assert abs(exp_alpha.energies - exp_alpha0.energies).max() < 1e-4

    # Check the hartree-fock energy
    dm = wfn.dm_alpha
    hf1 = sum([
        -2*hartree.expectation_value(dm),
        +1*exchange.expectation_value(dm),
    ]) + exp_alpha.energies[:wfn.nep].sum()*2
    hf2 = sum([
        2*cache['kin'].expectation_value(dm),
        -2*cache['na'].expectation_value(dm),
        +2*hartree.expectation_value(dm),
        -exchange.expectation_value(dm),
    ])
    enn = 9.2535672047 # nucleus-nucleus interaction
    assert abs(hf1 + enn - (-74.9592923284)) < 1e-4
    assert abs(hf2 + enn - (-74.9592923284)) < 1e-4


def test_dense_one_body_trace():
    lf = DenseLinalgFactory()
    op1 = lf.create_one_body(3)
    op1._array[:] = np.random.uniform(-1, 1, (3,3))
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


def test_dense_linalg_factory_properties():
    lf = DenseLinalgFactory(5)
    assert lf.default_nbasis == 5
    lf = DenseLinalgFactory()
    assert lf.default_nbasis is None
    lf.default_nbasis = 10
    ex = lf.create_expansion()
    assert ex.nbasis == 10
    assert ex.nfn == 10
    assert ex.energies.shape == (10,)
    assert ex.occupations.shape == (10,)
    op1 = lf.create_one_body()
    assert op1.nbasis == 10
    op2 = lf.create_two_body()
    assert op2.nbasis == 10


def test_dense_expansion_properties():
    lf = DenseLinalgFactory()
    ex = lf.create_expansion(10, 8)
    assert ex.nbasis == 10
    assert ex.nfn == 8
    assert ex.coeffs.shape == (10,8) # orbitals stored as columns
    assert ex.energies.shape == (8,)
    assert ex.occupations.shape == (8,)


def test_dense_one_body_properties():
    lf = DenseLinalgFactory()
    op = lf.create_one_body(3)
    assert op.nbasis == 3
    op.set_element(0, 1, 1.2)
    assert op.get_element(0, 1) == 1.2


def test_dense_two_body_properties():
    lf = DenseLinalgFactory()
    op = lf.create_two_body(3)
    assert op.nbasis == 3


def test_dense_one_body_assign():
    lf = DenseLinalgFactory()
    op1 = lf.create_one_body(3)
    op2 = lf.create_one_body(3)
    op1._array[:] = np.random.uniform(0, 1, (3, 3))
    op2.assign(op1)
    assert (op1._array == op2._array).all()


def test_dense_one_body_copy():
    lf = DenseLinalgFactory()
    op1 = lf.create_one_body(3)
    op1._array[:] = np.random.uniform(0, 1, (3, 3))
    op2 = op1.copy()
    assert (op1._array == op2._array).all()


def test_dense_expansion_copy():
    lf = DenseLinalgFactory()
    exp1 = lf.create_expansion(3, 2)
    exp1._coeffs[:] = np.random.uniform(0, 1, (3, 2))
    exp1._energies[:] = np.random.uniform(0, 1, 2)
    exp1._occupations[:] = np.random.uniform(0, 1, 2)
    exp2 = exp1.copy()
    assert (exp1._coeffs == exp2._coeffs).all()
    assert (exp1._energies == exp2._energies).all()
    assert (exp1._occupations == exp2._occupations).all()


def test_homo_lumo_ch3_hf():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    assert mol.wfn.exp_alpha.get_homo_index() == 4
    assert mol.wfn.exp_beta.get_homo_index() == 3
    assert mol.wfn.exp_alpha.get_lumo_index() == 5
    assert mol.wfn.exp_beta.get_lumo_index() == 4
    assert mol.wfn.exp_alpha.get_homo_index(1) == 3
    assert mol.wfn.exp_beta.get_homo_index(1) == 2
    assert mol.wfn.exp_alpha.get_lumo_index(1) == 6
    assert mol.wfn.exp_beta.get_lumo_index(1) == 5
    assert abs(mol.wfn.exp_alpha.get_homo_energy() - -3.63936540E-01) < 1e-8
    assert abs(mol.wfn.exp_alpha.get_homo_energy(1) - -5.37273275E-01) < 1e-8
    assert abs(mol.wfn.exp_alpha.get_lumo_energy() - 6.48361367E-01) < 1e-8
    assert abs(mol.wfn.exp_beta.get_homo_energy() - -5.18988806E-01) < 1e-8
    assert abs(mol.wfn.exp_beta.get_homo_energy(1) - -5.19454722E-01) < 1e-8
    assert abs(mol.wfn.exp_beta.get_lumo_energy() - 3.28562907E-01) < 1e-8
    assert abs(mol.wfn.exp_alpha.homo_energy - -3.63936540E-01) < 1e-8
    assert abs(mol.wfn.exp_alpha.lumo_energy - 6.48361367E-01) < 1e-8
    assert abs(mol.wfn.exp_beta.homo_energy - -5.18988806E-01) < 1e-8
    assert abs(mol.wfn.exp_beta.lumo_energy - 3.28562907E-01) < 1e-8


def test_naturals():
    fn_fchk = context.get_fn('test/ch3_hf_sto3g.fchk')
    mol = Molecule.from_file(fn_fchk)
    mol.wfn.clear_exp()
    exp_alpha = mol.wfn.init_exp('alpha')
    exp_alpha.derive_naturals(mol.wfn.dm_alpha, mol.obasis.compute_overlap(mol.lf))
    assert exp_alpha.occupations.min() > -1e-6
    assert exp_alpha.occupations.max() < 1+1e-6
    exp_alpha.check_normalization(mol.obasis.compute_overlap(mol.lf))
