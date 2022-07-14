# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2022 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


import numpy as np
from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def test_1d_hubbard_hamiltonian():
    # Test the half-filled 1-D Hubbard model Hamiltonian
    # with 10 sites for U=2 and using PBC

    occ_model = AufbauOccModel(5)
    modelham = Hubbard(10, pbc=True)
    orb_alpha = Orbitals(10)
    olp = modelham.compute_overlap()
    # t-param, t = -1
    kin = modelham.compute_kinetic(-1)
    # U-param, U = 2
    er = modelham.compute_er(2)
    # Guess
    guess_core_hamiltonian(olp, kin, orb_alpha)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
    ]
    ham = REffHam(terms)
    scf_solver = PlainSCFSolver()
    scf_solver(ham, olp, occ_model, orb_alpha)
    energy = ham.compute_energy()

    np.testing.assert_almost_equal(energy, -7.94427, decimal=4)


def test_1d_hubbard_hamiltonian_no_pbc():
    # Test the half-filled 1-D Hubbard model Hamiltonian
    # with 10 sites for U=2 without PBC

    occ_model = AufbauOccModel(5)
    modelham = Hubbard(10, pbc=False)
    orb_alpha = Orbitals(10)
    olp = modelham.compute_overlap()
    # t-param, t = -1
    kin = modelham.compute_kinetic(-1)
    # U-param, U = 2
    er = modelham.compute_er(2)
    # Guess
    guess_core_hamiltonian(olp, kin, orb_alpha)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
    ]
    ham = REffHam(terms)
    scf_solver = PlainSCFSolver()
    scf_solver(ham, olp, occ_model, orb_alpha)
    energy = ham.compute_energy()

    np.testing.assert_almost_equal(energy, -7.0533, decimal=4)
