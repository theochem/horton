# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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

import pytest

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_hf_cs_hf, check_lih_os_hf, \
    check_vanadium_sc_hf


def test_hf_cs_hf():
    check_hf_cs_hf(PlainSCFSolver(threshold=1e-10))


def test_lih_os_hf():
    check_lih_os_hf(PlainSCFSolver(threshold=1e-10))


def test_vanadium_sc_hf():
    with pytest.raises(NoSCFConvergence):
        check_vanadium_sc_hf(PlainSCFSolver(threshold=1e-10, maxiter=10))


def test_hf_cs_hf_level_shift():
    check_hf_cs_hf(PlainSCFSolver(threshold=1e-10, level_shift=0.01))


def test_lih_os_hf_level_shift():
    check_lih_os_hf(PlainSCFSolver(threshold=1e-10, level_shift=0.02))


def test_vanadium_sc_hf_level_shift():
    with pytest.raises(ValueError):
        check_vanadium_sc_hf(PlainSCFSolver(threshold=1e-10, level_shift=-0.1))


def test_hf_water_321g_mistake():
    # When one forgets to construct the initial guess, some error must be
    # raised...
    fn_xyz = context.get_fn('test/water.xyz')
    mol = IOData.from_file(fn_xyz)
    obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')
    occ_model = AufbauOccModel(5)
    orb_alpha = Orbitals(obasis.nbasis)
    olp = obasis.compute_overlap()
    kin = obasis.compute_kinetic()
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers)
    er = obasis.compute_electron_repulsion()
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)
    scf_solver = PlainSCFSolver()
    with pytest.raises(AssertionError):
        scf_solver(ham, olp, occ_model, orb_alpha)
