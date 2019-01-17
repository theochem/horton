# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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

from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_hf_cs_hf, check_lih_os_hf, \
    check_vanadium_sc_hf


def test_hf_cs_hf():
    check_hf_cs_hf(PlainSCFSolver(threshold=1e-10))


def test_lih_os_hf():
    check_lih_os_hf(PlainSCFSolver(threshold=1e-10))


def test_vanadium_sc_hf():
    with assert_raises(NoSCFConvergence):
        check_vanadium_sc_hf(PlainSCFSolver(threshold=1e-10, maxiter=10))


def test_hf_cs_hf_level_shift():
    check_hf_cs_hf(PlainSCFSolver(threshold=1e-10, level_shift=0.01))


def test_lih_os_hf_level_shift():
    check_lih_os_hf(PlainSCFSolver(threshold=1e-10, level_shift=0.02))


def test_vanadium_sc_hf_level_shift():
    with assert_raises(ValueError):
        check_vanadium_sc_hf(PlainSCFSolver(threshold=1e-10, level_shift=-0.1))


def test_hf_water_321g_mistake():
    # When one forgets to construct the initial guess, some error must be
    # raised...
    fn_xyz = context.get_fn('test/water.xyz')
    mol = IOData.from_file(fn_xyz)
    obasis = get_gobasis(mol.coordinates, mol.numbers, '3-21G')
    lf = DenseLinalgFactory(obasis.nbasis)
    occ_model = AufbauOccModel(5)
    exp_alpha = lf.create_expansion(obasis.nbasis)
    olp = obasis.compute_overlap(lf)
    kin = obasis.compute_kinetic(lf)
    na = obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, lf)
    er = obasis.compute_electron_repulsion(lf)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)
    scf_solver = PlainSCFSolver()
    with assert_raises(AssertionError):
        scf_solver(ham, lf, olp, occ_model, exp_alpha)
