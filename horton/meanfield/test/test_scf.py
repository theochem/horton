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
    sys = System.from_file(fn_fchk)

    olp = sys.get_overlap()
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    external = {'nn': compute_nucnuc(sys.coordinates, sys.numbers)}
    terms = [
        OneBodyTerm(kin, sys.wfn, 'kin'),
        DirectTerm(er, sys.wfn),
        ExchangeTerm(er, sys.wfn),
        OneBodyTerm(nai, sys.wfn, 'ne'),
    ]
    ham = Hamiltonian(terms, external)

    guess_core_hamiltonian(sys.wfn, olp, kin, nai)
    assert convergence_error_eigen(ham, sys.wfn, sys.lf, olp) > 1e-8
    converge_scf(ham, sys.wfn, sys.lf, olp)
    assert convergence_error_eigen(ham, sys.wfn, sys.lf, olp) < 1e-8

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
    assert abs(sys.wfn.exp_alpha.energies - expected_alpha_energies).max() < 1e-5
    assert abs(sys.wfn.exp_beta.energies - expected_beta_energies).max() < 1e-5

    ham.compute()
    # compare with g09
    assert abs(ham.cache['energy'] - -7.687331212191962E+00) < 1e-8
    assert abs(ham.cache['energy_kin'] - 7.640603924034E+00) < 2e-7
    assert abs(ham.cache['energy_hartree'] + ham.cache['energy_exchange_hartree_fock'] - 2.114420907894E+00) < 1e-7
    assert abs(ham.cache['energy_ne'] - -1.811548789281E+01) < 2e-7
    assert abs(ham.cache['energy_nn'] - 0.6731318487) < 1e-8


def test_hf_water_321g_mistake():
    fn_xyz = context.get_fn('test/water.xyz')
    sys = System.from_file(fn_xyz, obasis='3-21G')
    wfn = setup_mean_field_wfn(sys.obasis.nbasis, sys.pseudo_numbers, sys.lf, charge=0)
    olp = sys.get_overlap()
    kin = sys.get_kinetic()
    nai = sys.get_nuclear_attraction()
    er = sys.get_electron_repulsion()
    terms = [
        OneBodyTerm(kin, wfn, 'kin'),
        DirectTerm(er, wfn),
        ExchangeTerm(er, wfn),
        OneBodyTerm(nai, wfn, 'ne'),
    ]
    ham = Hamiltonian(terms)
    with assert_raises(AttributeError):
        converge_scf(ham, wfn, sys.lf, olp)
