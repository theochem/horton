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


def test_fock_n2_hfs_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    grid = get_some_grid(sys)
    libxc_term = LibXCLDA('x')
    ham1 = Hamiltonian(sys, [Hartree(), libxc_term], grid)
    builtin_term = DiracExchange()
    ham2 = Hamiltonian(sys, [Hartree(), builtin_term], grid)

    # Compare the potential computed by libxc with the builtin implementation
    libxc_term._update_operator()
    libxc_pot = libxc_term.cache.load('pot_libxc_lda_x_alpha')
    builtin_term._update_exchange()
    builtin_pot = builtin_term.cache.load('pot_exchange_dirac_alpha')
    rho = libxc_term.update_rho('alpha')
    # Libxc apparently approximates values of the potential below 1e-4 with zero.
    assert abs(libxc_pot - builtin_pot).max() < 1e-4

    # Check of the libxc energy matches our implementation
    energy1 = ham1.compute_energy()
    ex1 = sys.props['energy_libxc_lda_x']
    energy2 = ham2.compute_energy()
    ex2 = sys.props['energy_exchange_dirac']
    assert abs(ex1 - ex2) < 1e-10
    assert abs(energy1 - energy2) < 1e-10

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error(ham1) < 1e-5

    # Converge from scratch
    guess_hamiltonian_core(sys)
    assert convergence_error(ham1) > 1e-8
    assert converge_scf(ham1)
    assert convergence_error(ham1) < 1e-8

    # test orbital energies
    expected_energies = np.array([
        -1.37107053E+01, -1.37098006E+01, -9.60673085E-01, -3.57928483E-01,
        -3.16017655E-01, -3.16017655E-01, -2.12998316E-01, 6.84030479E-02,
        6.84030479E-02, 7.50192517E-01,
    ])
    assert abs(sys.wfn.exp_alpha.energies - expected_energies).max() < 3e-5

    ham1.compute_energy()
    # compare with g09
    assert abs(sys.props['energy_ne'] - -2.981579553570E+02) < 1e-5
    assert abs(sys.props['energy_kin'] - 1.061620887711E+02) < 1e-5
    assert abs(sys.props['energy_hartree'] + sys.props['energy_libxc_lda_x'] - 6.247259253877E+01) < 1e-4
    assert abs(sys.props['energy'] - -106.205213597) < 1e-4
    assert abs(sys.props['energy_nn'] - 23.3180604505) < 1e-8


def test_hamiltonian_h3_hfs_321g():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)
    grid = get_some_grid(sys)
    libxc_term = LibXCLDA('x')
    ham1 = Hamiltonian(sys, [Hartree(), libxc_term], grid)
    builtin_term = DiracExchange()
    ham2 = Hamiltonian(sys, [Hartree(), builtin_term], grid)

    # Compare the potential computed by libxc with the builtin implementation
    libxc_term._update_operator()
    libxc_pot = libxc_term.cache.load('pot_libxc_lda_x_both')[:,0]
    builtin_term._update_exchange()
    builtin_pot = builtin_term.cache.load('pot_exchange_dirac_alpha')
    rho = libxc_term.update_rho('alpha')
    # Libxc apparently approximates values of the potential below 1e-4 with zero.
    assert abs(libxc_pot - builtin_pot).max() < 1e-4

    # Check of the libxc energy matches our implementation
    energy1 = ham1.compute_energy()
    ex1 = sys.props['energy_libxc_lda_x']
    energy2 = ham2.compute_energy()
    ex2 = sys.props['energy_exchange_dirac']
    assert abs(ex1 - ex2) < 1e-10
    assert abs(energy1 - energy2) < 1e-10

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error(ham1) < 1e-5

    # Converge from scratch
    guess_hamiltonian_core(sys)
    assert convergence_error(ham1) > 1e-8
    assert converge_scf_oda(ham1)
    assert convergence_error(ham1) < 1e-8

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

    ham1.compute_energy()
    # compare with g09
    assert abs(sys.props['energy_ne'] - -6.832069993374E+00) < 1e-5
    assert abs(sys.props['energy_kin'] - 1.870784279014E+00) < 1e-5
    assert abs(sys.props['energy_hartree'] + sys.props['energy_libxc_lda_x'] - 1.658810998195E+00) < 1e-6
    assert abs(sys.props['energy'] - -1.412556114057104E+00) < 1e-5
    assert abs(sys.props['energy_nn'] - 1.8899186021) < 1e-8


def test_co_pbe_sto3g():
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    grid = get_some_grid(sys)
    libxc_x_term = LibXCGGA('x_pbe')
    libxc_c_term = LibXCGGA('c_pbe')
    ham = Hamiltonian(sys, [Hartree(), libxc_x_term, libxc_c_term], grid)

    # Test energy before scf
    ham.compute_energy()
    assert abs(sys.props['energy'] - -1.116465967841901E+02) < 1e-4

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error(ham) < 1e-5

    # Converge from scratch
    guess_hamiltonian_core(sys)
    assert convergence_error(ham) > 1e-5
    assert converge_scf_oda(ham, threshold=1e-3)
    assert convergence_error(ham) < 1e-5

    # test orbital energies
    expected_energies = np.array([
         -1.86831122E+01, -9.73586915E+00, -1.03946082E+00, -4.09331776E-01,
         -3.48686522E-01, -3.48686522E-01, -2.06049056E-01, 5.23730418E-02,
         5.23730418E-02, 6.61093726E-01
    ])
    assert abs(sys.wfn.exp_alpha.energies - expected_energies).max() < 1e-2

    ham.compute_energy()
    # compare with g09
    assert abs(sys.props['energy_ne'] - -3.072370116827E+02) < 1e-2
    assert abs(sys.props['energy_kin'] - 1.103410779827E+02) < 1e-2
    assert abs(sys.props['energy_hartree'] + sys.props['energy_libxc_gga_x_pbe'] + sys.props['energy_libxc_gga_c_pbe'] - 6.273115782683E+01) < 1e-2
    assert abs(sys.props['energy'] - -1.116465967841901E+02) < 1e-4
    assert abs(sys.props['energy_nn'] - 22.5181790889) < 1e-7


def test_h3_pbe_321g():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    sys = System.from_file(fn_fchk)
    grid = get_some_grid(sys)
    libxc_x_term = LibXCGGA('x_pbe')
    libxc_c_term = LibXCGGA('c_pbe')
    ham = Hamiltonian(sys, [Hartree(), libxc_x_term, libxc_c_term], grid)

    # compute the energy before converging
    ham.compute_energy()
    assert abs(sys.props['energy'] - -1.593208400939354E+00) < 1e-5

    # The convergence should be reasonable, not perfect because of limited
    # precision in Gaussian fchk file:
    assert convergence_error(ham) < 2e-6

    # Converge from scratch
    guess_hamiltonian_core(sys)
    assert convergence_error(ham) > 1e-5
    assert converge_scf_oda(ham, threshold=1e-5)
    assert convergence_error(ham) < 1e-5

    # test orbital energies
    expected_energies = np.array([
        -5.41141676E-01, -1.56826691E-01, 2.13089637E-01, 7.13565167E-01,
        7.86810564E-01, 1.40663544E+00
    ])
    assert abs(sys.wfn.exp_alpha.energies - expected_energies).max() < 2e-5
    expected_energies = np.array([
        -4.96730336E-01, -5.81411249E-02, 2.73586652E-01, 7.41987185E-01,
        8.76161160E-01, 1.47488421E+00
    ])
    assert abs(sys.wfn.exp_beta.energies - expected_energies).max() < 2e-5

    ham.compute_energy()
    # compare with g09
    assert abs(sys.props['energy_ne'] - -6.934705182067E+00) < 1e-5
    assert abs(sys.props['energy_kin'] - 1.948808793424E+00) < 1e-5
    assert abs(sys.props['energy_hartree'] + sys.props['energy_libxc_gga_x_pbe'] + sys.props['energy_libxc_gga_c_pbe'] - 1.502769385597E+00) < 1e-5
    assert abs(sys.props['energy'] - -1.593208400939354E+00) < 1e-5
    assert abs(sys.props['energy_nn'] - 1.8899186021) < 1e-8


def test_cubic_interpolation_c_pbe_cs():
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), LibXCGGA('c_pbe')], grid)

    dm0 = sys.wfn.dm_alpha.copy()
    converge_scf_oda(ham, max_iter=1)
    dm1 = sys.wfn.dm_alpha.copy()

    check_cubic_cs_wrapper(ham, dm0, dm1)


def test_cubic_interpolation_x_pbe_cs():
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), LibXCGGA('x_pbe')], grid)

    dm0 = sys.wfn.dm_alpha.copy()
    converge_scf_oda(ham, max_iter=1)
    dm1 = sys.wfn.dm_alpha.copy()

    check_cubic_cs_wrapper(ham, dm0, dm1)


def test_cubic_interpolation_hfs_cs():
    fn_fchk = context.get_fn('test/co_pbe_sto3g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), LibXCLDA('x')], grid)

    dm0 = sys.wfn.dm_alpha.copy()
    converge_scf_oda(ham, max_iter=1)
    dm1 = sys.wfn.dm_alpha.copy()

    check_cubic_cs_wrapper(ham, dm0, dm1)


def test_cubic_interpolation_o3lyp_cs():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    libxc_term = LibXCHybridGGA('xc_o3lyp')
    ham = Hamiltonian(sys, [HartreeFockExchange(libxc_term.get_exx_fraction()), libxc_term], grid)

    dm0 = sys.wfn.dm_alpha.copy()
    converge_scf_oda(ham, max_iter=1)
    dm1 = sys.wfn.dm_alpha.copy()

    check_cubic_cs_wrapper(ham, dm0, dm1)


def test_cubic_interpolation_c_pbe_os():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), LibXCGGA('c_pbe')], grid)

    dma0 = sys.wfn.dm_alpha.copy()
    dmb0 = sys.wfn.dm_beta.copy()
    converge_scf_oda(ham, max_iter=1)
    dma1 = sys.wfn.dm_alpha.copy()
    dmb1 = sys.wfn.dm_beta.copy()

    check_cubic_os_wrapper(ham, dma0, dmb0, dma1, dmb1)


def test_cubic_interpolation_x_pbe_os():
    fn_fchk = context.get_fn('test/h3_pbe_321g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), LibXCGGA('x_pbe')], grid)

    dma0 = sys.wfn.dm_alpha.copy()
    dmb0 = sys.wfn.dm_beta.copy()
    converge_scf_oda(ham, max_iter=1)
    dma1 = sys.wfn.dm_alpha.copy()
    dmb1 = sys.wfn.dm_beta.copy()

    check_cubic_os_wrapper(ham, dma0, dmb0, dma1, dmb1)


def test_cubic_interpolation_hfs_os():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    ham = Hamiltonian(sys, [Hartree(), LibXCLDA('x')], grid)

    dma0 = sys.wfn.dm_alpha.copy()
    dmb0 = sys.wfn.dm_beta.copy()
    guess_hamiltonian_core(sys)
    dma1 = sys.wfn.dm_alpha.copy()
    dmb1 = sys.wfn.dm_beta.copy()

    check_cubic_os_wrapper(ham, dma0, dmb0, dma1, dmb1)


def test_cubic_interpolation_o3lyp_os():
    fn_fchk = context.get_fn('test/h3_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)

    grid = get_some_grid(sys)
    libxc_term = LibXCHybridGGA('xc_o3lyp')
    ham = Hamiltonian(sys, [HartreeFockExchange(libxc_term.get_exx_fraction()), libxc_term], grid)

    dma0 = sys.wfn.dm_alpha.copy()
    dmb0 = sys.wfn.dm_beta.copy()
    guess_hamiltonian_core(sys)
    dma1 = sys.wfn.dm_alpha.copy()
    dmb1 = sys.wfn.dm_beta.copy()

    check_cubic_os_wrapper(ham, dma0, dmb0, dma1, dmb1)


def test_hyb_gga_exx_fraction():
    t = LibXCHybridGGA('xc_pbeh') # The PBE0 functional
    assert t.get_exx_fraction() == 0.25


def test_lda_c_vwn_present():
    t = LibXCLDA('c_vwn')     # The VWN 5 functional
    t = LibXCLDA('c_vwn_4')   # The VWN 4 functional
