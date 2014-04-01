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
from horton.meanfield.test.common import check_scf_hf_cs_hf


def test_scf_cs_hf():
    check_scf_hf_cs_hf(SCFWrapper('plain', threshold=1e-10))


def test_scf_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = Molecule.from_file(fn_fchk)

    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    nai = mol.obasis.compute_nuclear_attraction(mol.pseudo_numbers, mol.coordinates, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    external = {'nn': compute_nucnuc(mol.coordinates, mol.numbers)}
    terms = [
        OneBodyTerm(kin, mol.wfn, 'kin'),
        DirectTerm(er, mol.wfn),
        ExchangeTerm(er, mol.wfn),
        OneBodyTerm(nai, mol.wfn, 'ne'),
    ]
    ham = Hamiltonian(terms, external)

    guess_core_hamiltonian(mol.wfn, olp, kin, nai)
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) > 1e-8
    converge_scf(ham, mol.wfn, mol.lf, olp)
    assert convergence_error_eigen(ham, mol.wfn, mol.lf, olp) < 1e-8

    expected_alpha_energies = np.array([
        -2.76116635E+00, -7.24564188E-01, -1.79148636E-01, -1.28235698E-01,
        -1.28235698E-01, -7.59817520E-02, -1.13855167E-02, 6.52484445E-03,
        6.52484445E-03, 7.52201895E-03, 9.70893294E-01,
    ])
    expected_beta_energies = np.array([
        -2.76031162E+00, -2.08814026E-01, -1.53071066E-01, -1.25264964E-01,
        -1.25264964E-01, -1.24605870E-02, 5.12761388E-03, 7.70499854E-03,
        7.70499854E-03, 2.85176080E-02, 1.13197479E+00,
    ])
    assert abs(mol.wfn.exp_alpha.energies - expected_alpha_energies).max() < 1e-5
    assert abs(mol.wfn.exp_beta.energies - expected_beta_energies).max() < 1e-5

    ham.compute()
    # compare with g09
    assert abs(ham.cache['energy'] - -7.687331212191962E+00) < 1e-8
    assert abs(ham.cache['energy_kin'] - 7.640603924034E+00) < 2e-7
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_exchange_hartree_fock'] - 2.114420907894E+00) < 1e-7
    assert abs(ham.cache['energy_ne'] - -1.811548789281E+01) < 2e-7
    assert abs(ham.cache['energy_nn'] - 0.6731318487) < 1e-8


def test_hf_water_321g_mistake():
    fn_xyz = context.get_fn('test/water.xyz')
    mol = Molecule.from_file(fn_xyz)
    obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')
    wfn = setup_mean_field_wfn(obasis.nbasis, mol.numbers, mol.lf, charge=0)
    lf = DenseLinalgFactory(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    nai = obasis.compute_nuclear_attraction(mol.pseudo_numbers, mol.coordinates, lf)
    er = obasis.compute_electron_repulsion(lf)
    terms = [
        OneBodyTerm(kin, wfn, 'kin'),
        DirectTerm(er, wfn),
        ExchangeTerm(er, wfn),
        OneBodyTerm(nai, wfn, 'ne'),
    ]
    ham = Hamiltonian(terms)
    with assert_raises(AttributeError):
        converge_scf(ham, wfn, lf, olp)
