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


from nose.tools import assert_raises

from .common import check_hf_cs_hf, check_lih_os_hf, \
    load_olp, load_kin, load_na, load_er
from .. import PlainSCFSolver, AufbauOccModel, Orbitals, RTwoIndexTerm, \
    RDirectTerm, RExchangeTerm, REffHam


def test_hf_cs_hf():
    check_hf_cs_hf(PlainSCFSolver(threshold=1e-10))


def test_lih_os_hf():
    check_lih_os_hf(PlainSCFSolver(threshold=1e-10))


# TODO: Move to higher level test
# def test_vanadium_sc_hf():
#     with assert_raises(NoSCFConvergence):
#         check_vanadium_sc_hf(PlainSCFSolver(threshold=1e-10, maxiter=10))


def test_hf_cs_hf_level_shift():
    check_hf_cs_hf(PlainSCFSolver(threshold=1e-10, level_shift=0.01))


def test_lih_os_hf_level_shift():
    check_lih_os_hf(PlainSCFSolver(threshold=1e-10, level_shift=0.02))


# TODO: Move to higher level test
# def test_vanadium_sc_hf_level_shift():
#     with assert_raises(ValueError):
#         check_vanadium_sc_hf(PlainSCFSolver(threshold=1e-10, level_shift=-0.1))


def test_hf_water_321g_mistake():
    # When one forgets to construct the initial guess, some error must be
    # raised...
    fname = 'water_3_21G_xyz'
    occ_model = AufbauOccModel(5)
    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    terms = [
        RTwoIndexTerm(kin, 'kin'),
        RDirectTerm(er, 'hartree'),
        RExchangeTerm(er, 'x_hf'),
        RTwoIndexTerm(na, 'ne'),
    ]
    ham = REffHam(terms)
    scf_solver = PlainSCFSolver()
    orb_alpha = Orbitals(olp.shape[0])
    with assert_raises(AssertionError):
        scf_solver(ham, olp, occ_model, orb_alpha)
