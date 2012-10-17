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


def test_fock_n2_hfs_sto3g():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    int1d = TrapezoidIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    grid = BeckeMolGrid(sys, (rtf, int1d, 110), random_rotate=False)
    libxc_term = LibXCTerm('lda_x')
    ham1 = Hamiltonian(sys, [Hartree(), libxc_term], grid)
    builtin_term = DiracExchange()
    ham2 = Hamiltonian(sys, [Hartree(), builtin_term], grid)

    # Compare the potential computed by libxc with the builtin implementation
    libxc_term._update_operator()
    libxc_pot = libxc_term.cache.load('pot_libxc_lda_x_alpha')
    builtin_term._update_exchange()
    builtin_pot = builtin_term.cache.load('pot_exchange_dirac_alpha')
    rho = libxc_term.update_rho('alpha')
    print libxc_pot - builtin_pot
    print abs(libxc_pot - builtin_pot).max()
    # Libxc apparently approximates values of the potential below 1e-4 with zero.
    assert abs(libxc_pot - builtin_pot).max() < 1e-4

    # Check of the libxc energy matches our implementation
    energy1 = ham1.compute_energy()
    ex1 = sys.props['energy_libxc_lda_x']
    energy2 = ham2.compute_energy()
    ex2 = sys.props['energy_exchange_dirac']
    print ex1, ex2
    assert abs(ex1 - ex2) < 1e-10
    print energy1, energy2
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
    assert abs(sys.props['energy_hartree'] + sys.props['energy_exchange_dirac'] - 6.247259253877E+01) < 1e-4
    assert abs(sys.props['energy'] - -106.205213597) < 1e-4
    assert abs(sys.props['energy_nn'] - 23.3180604505) < 1e-8
