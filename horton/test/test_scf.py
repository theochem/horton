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

    guess_hamiltonian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    assert convergence_error(ham) > 1e-8
    assert converge_scf(ham)
    assert convergence_error(ham) < 1e-8

    # get the hatreefock term
    for term in ham.terms:
        if isinstance(term, HartreeFock):
            hf_term = term
            break

    # test operator consistency
    my_coulomb = sys.lf.create_one_body()
    dm_alpha = sys.wfn.dm_alpha
    hf_term.electron_repulsion.apply_direct(dm_alpha, my_coulomb)
    my_coulomb.iscale(2)
    error = abs(my_coulomb._array - ham.cache.load('op_coulomb')._array).max()
    assert error < 1e-5

    # test orbital energies
    expected_energies = np.array([
        -2.59083334E+01, -1.44689996E+00, -5.57467136E-01, -4.62288194E-01,
        -4.62288194E-01, 5.39578910E-01,
    ])
    assert abs(sys.wfn.exp_alpha.energies - expected_energies).max() < 1e-5

    ham.compute_energy()
    # compare with g09
    assert abs(sys.props['energy'] - -9.856961609951867E+01) < 1e-8
    assert abs(sys.props['energy_kin'] - 9.766140786239E+01) < 2e-7
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_fock'] - 4.561984106482E+01) < 1e-7
    assert abs(sys.props['energy_ne'] - -2.465756615329E+02) < 2e-7
    assert abs(sys.props['energy_nn'] - 4.7247965053) < 1e-8


def test_scf_os():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)

    guess_hamiltonian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    assert convergence_error(ham) > 1e-8
    assert converge_scf(ham)
    assert convergence_error(ham) < 1e-8

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

    ham.compute_energy()
    # compare with g09
    assert abs(sys.props['energy'] - -7.687331212191962E+00) < 1e-8
    assert abs(sys.props['energy_kin'] - 7.640603924034E+00) < 2e-7
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_fock'] - 2.114420907894E+00) < 1e-7
    assert abs(sys.props['energy_ne'] - -1.811548789281E+01) < 2e-7
    assert abs(sys.props['energy_nn'] - 0.6731318487) < 1e-8


def test_scf_oda_water_hfs_321g():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)

    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 2e1, 100)
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False)
    ham = Hamiltonian(sys, [Hartree(), DiracExchange()], grid)

    if True:
        # The convergence should be reasonable, not perfect because of limited
        # precision in Gaussian fchk file:
        assert convergence_error(ham) < 2e-5

        # The energies should also be in reasonable agreement. Repeated to check for
        # stupid bugs
        for i in xrange(2):
            ham.invalidate()
            ham.compute_energy()
            expected_energies = np.array([
                -1.83691041E+01, -8.29412411E-01, -4.04495188E-01, -1.91740814E-01,
                -1.32190590E-01, 1.16030419E-01, 2.08119657E-01, 9.69825207E-01,
                9.99248500E-01, 1.41697384E+00, 1.47918828E+00, 1.61926596E+00,
                2.71995350E+00
            ])

            assert abs(sys.wfn.exp_alpha.energies - expected_energies).max() < 2e-4
            assert abs(sys.props['energy_ne'] - -1.977921986200E+02) < 1e-7
            assert abs(sys.props['energy_kin'] - 7.525067610865E+01) < 1e-9
            assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_dirac'] - 3.864299848058E+01) < 1e-4
            assert abs(sys.props['energy'] - -7.474134898935590E+01) < 1e-4
            assert abs(sys.props['energy_nn'] - 9.1571750414) < 2e-8

    # Converge from scratch
    guess_hamiltonian_core(sys)
    assert convergence_error(ham) > 1e-5
    assert converge_scf_oda(ham, threshold=1e-5)
    assert convergence_error(ham) < 1e-5

    assert abs(sys.props['energy_ne'] - -1.977921986200E+02) < 1e-4
    assert abs(sys.props['energy_kin'] - 7.525067610865E+01) < 3e-5
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_dirac'] - 3.864299848058E+01) < 1e-4
    assert abs(sys.props['energy'] - -7.474134898935590E+01) < 1e-4


def test_scf_oda_water_hf_321g():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)
    ham = Hamiltonian(sys, [HartreeFock()])

    # test continuation of interupted scf_oda
    guess_hamiltonian_core(sys)
    e0 = ham.compute_energy()
    assert convergence_error(ham) > 1e-5
    assert not converge_scf_oda(ham, threshold=1e-2, max_iter=3)
    assert not sys.wfn._cache.has('exp_alpha')
    e1 = ham.compute_energy()
    assert not converge_scf_oda(ham, threshold=1e-2, max_iter=3)
    e2 = ham.compute_energy()
    assert e1 < e0
    assert e2 < e1


def test_scf_oda_lih_hfs_321g():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 2e1, 100)
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False)
    ham = Hamiltonian(sys, [Hartree(), DiracExchange()], grid)

    # test continuation of interupted scf_oda
    guess_hamiltonian_core(sys)
    e0 = ham.compute_energy()
    assert convergence_error(ham) > 1e-5
    assert not converge_scf_oda(ham, threshold=1e-3, max_iter=3)
    assert not sys.wfn._cache.has('exp_alpha')
    e1 = ham.compute_energy()
    assert not converge_scf_oda(ham, threshold=1e-3, max_iter=3)
    e2 = ham.compute_energy()
    assert e1 < e0
    assert e2 < e1

    # continue till convergence
    assert converge_scf_oda(ham, threshold=1e-3)


def test_hf_water_321g_mistake():
    fn_xyz = context.get_fn('test/water.xyz')
    sys = System.from_file(fn_xyz, obasis='3-21G')
    sys.init_wfn(charge=0)
    ham = Hamiltonian(sys, [HartreeFock()])
    try:
        converge_scf(ham)
        assert False
    except AttributeError:
        pass


def test_find_min_cubic():
    from horton.scf import find_min_cubic
    assert find_min_cubic(0.2, 0.5, 3.0, -0.7) == 0.0
    assert abs(find_min_cubic(2.1, -5.2, -3.0, 2.8) - 0.939645667705) < 1e-8
    assert abs(find_min_cubic(0.0, 1.0, -0.1, -0.1) - 0.0153883154024) < 1e-8
    assert find_min_cubic(1.0, 1.0, 1.0, -1.0) == 1.0
    assert find_min_cubic(1.0, 1.0, -1.0, 1.0) == 0.5
    assert find_min_cubic(0.0, 1.0, 1.0, 1.0) == 0.0
    assert find_min_cubic(1.0, 0.0, -1.0, -1.0) == 1.0
    assert find_min_cubic(0.0, 1.0, 0.0, 0.0) == 0.0
    assert find_min_cubic(0.0, 1.0, 0.1, 0.1) == 0.0


def test_scf_oda_aufbau_spin():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.occ_model = AufbauSpinOccModel(3)

    guess_hamiltonian_core(sys)
    ham = Hamiltonian(sys, [HartreeFock()])
    assert converge_scf_oda(ham)
