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


def test_scf_cs():
    fn_fchk = context.get_fn('test/hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    guess_hamiltionian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    assert converge_scf(ham)

    # get the hatreefock term
    for term in ham.terms:
        if isinstance(term, HartreeFock):
            hf_term = term
            break

    # test operator consistency
    coulomb = sys.lf.create_one_body(sys.obasis.nbasis)
    hf_term.electron_repulsion.apply_direct(sys.dms['alpha'], coulomb)
    coulomb.iscale(2)
    error = abs(coulomb._array - hf_term.coulomb._array).max()
    assert error < 1e-5

    # test orbital energies
    expected_energies = np.array([
        -2.59083334E+01, -1.44689996E+00, -5.57467136E-01, -4.62288194E-01,
        -4.62288194E-01, 5.39578910E-01,
    ])
    assert abs(sys.wfn.expansion.energies - expected_energies).max() < 1e-5

    ham.compute_energy()
    # compare with g09
    assert abs(sys.props['energy'] - -9.856961609951867E+01) < 1e-8
    assert abs(sys.props['energy_kin'] - 9.766140786239E+01) < 2e-7
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_fock'] - 4.561984106482E+01) < 1e-7
    assert abs(sys.props['energy_ne'] - -2.465756615329E+02) < 2e-7
    assert abs(sys.props['energy_nn'] - 4.7247965053) < 1e-8

    # ugly hack:
    hf_term.exchange_beta = hf_term.exchange_alpha
    sys.dms['beta'] = sys.dms['alpha']
    dm_full = sys.lf.create_one_body(sys.obasis.nbasis)
    dm_full.iadd(sys.dms['alpha'], factor=2)
    sys.dms['full'] = dm_full
    ham.compute_energy()
    assert abs(sys.props['energy'] - -9.856961609951867E+01) < 1e-8
    assert abs(sys.props['energy_kin'] - 9.766140786239E+01) < 2e-7
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_fock'] - 4.561984106482E+01) < 1e-7
    assert abs(sys.props['energy_ne'] - -2.465756615329E+02) < 2e-7
    assert abs(sys.props['energy_nn'] - 4.7247965053) < 1e-8



def test_scf_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    guess_hamiltionian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    assert converge_scf(ham)
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
    assert abs(sys.wfn.alpha_expansion.energies - expected_alpha_energies).max() < 1e-5
    assert abs(sys.wfn.beta_expansion.energies - expected_beta_energies).max() < 1e-5

    ham.compute_energy()
    # compare with g09
    assert abs(sys.props['energy'] - -7.687331212191962E+00) < 1e-8
    assert abs(sys.props['energy_kin'] - 7.640603924034E+00) < 2e-7
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_fock'] - 2.114420907894E+00) < 1e-7
    assert abs(sys.props['energy_ne'] - -1.811548789281E+01) < 2e-7
    assert abs(sys.props['energy_nn'] - 0.6731318487) < 1e-8
