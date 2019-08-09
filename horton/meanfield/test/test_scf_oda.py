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
from nose.plugins.attrib import attr
from nose.tools import assert_raises

from .common import check_hf_cs_hf, check_lih_os_hf, \
    check_water_cs_hfs, check_n2_cs_hfs, check_h3_os_hfs, check_h3_os_pbe, \
    check_co_cs_pbe, check_water_cs_m05, \
    check_methyl_os_tpss, load_mdata, load_olp, load_kin, load_na, load_er, load_orbs_alpha, \
    load_orbs_beta
from .. import ODASCFSolver, AufbauSpinOccModel, UTwoIndexTerm, UDirectTerm, \
    UExchangeTerm, UEffHam, guess_core_hamiltonian, check_dm


def test_hf_cs_hf():
    check_hf_cs_hf(ODASCFSolver(threshold=1e-7))


def test_lih_os_hf():
    check_lih_os_hf(ODASCFSolver(threshold=1e-7))


def test_water_cs_hfs():
    check_water_cs_hfs(ODASCFSolver(threshold=1e-6))


@attr('slow')
def test_n2_cs_hfs():
    check_n2_cs_hfs(ODASCFSolver(threshold=1e-6))


@attr('slow')
def test_h3_os_hfs():
    check_h3_os_hfs(ODASCFSolver(threshold=1e-6))


@attr('slow')
def test_co_cs_pbe():
    check_co_cs_pbe(ODASCFSolver(threshold=1e-5))


@attr('slow')
def test_h3_os_pbe():
    check_h3_os_pbe(ODASCFSolver(threshold=1e-6))


# TODO: Move to higher level test
# def test_vanadium_sc_hf():
#     with assert_raises(NoSCFConvergence):
#         check_vanadium_sc_hf(ODASCFSolver(threshold=1e-10, maxiter=10))


def test_water_cs_m05():
    check_water_cs_m05(ODASCFSolver(threshold=1e-6))


def test_methyl_os_tpss():
    check_methyl_os_tpss(ODASCFSolver(threshold=1e-4))


def test_find_min_cubic():
    from ..scf_oda import find_min_cubic
    assert find_min_cubic(0.2, 0.5, 3.0, -0.7) == 0.0
    assert abs(find_min_cubic(2.1, -5.2, -3.0, 2.8) - 0.939645667705) < 1e-8
    assert abs(find_min_cubic(0.0, 1.0, -0.1, -0.1) - 0.0153883154024) < 1e-8
    assert find_min_cubic(1.0, 1.0, 1.0, -1.0) == 1.0
    assert find_min_cubic(1.0, 1.0, -1.0, 1.0) == 0.5
    assert find_min_cubic(0.0, 1.0, 1.0, 1.0) == 0.0
    assert find_min_cubic(1.0, 0.0, -1.0, -1.0) == 1.0
    assert find_min_cubic(0.0, 1.0, 0.0, 0.0) == 0.0
    assert find_min_cubic(0.0, 1.0, 0.1, 0.1) == 0.0


def test_find_min_quadratic():
    from ..scf_oda import find_min_quadratic
    assert find_min_quadratic(0.0, -0.7) == 1.0
    assert abs(find_min_quadratic(-3.0, 2.8) - 0.51724137931) < 1e-8
    assert abs(find_min_quadratic(-0.2, 0.1) - 0.666666666667) < 1e-8
    assert find_min_quadratic(-1.0, 1.0) == 0.5
    assert find_min_quadratic(-1.0, -1.0) == 1.0
    assert find_min_quadratic(1.0, 1.0) == 0.0
    assert find_min_quadratic(1.0, -0.5) == 0.0
    assert find_min_quadratic(1.0, -1.5) == 1.0


def test_aufbau_spin():
    fname = 'li_h_3_21G_hf_g09_fchk'
    mdata = load_mdata(fname)
    occ_model = AufbauSpinOccModel(3)

    olp = load_olp(fname)
    kin = load_kin(fname)
    na = load_na(fname)
    er = load_er(fname)
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er, 'x_hf'),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham = UEffHam(terms)

    # Construct an initial state with unstable spin polarization
    orb_alpha = load_orbs_alpha(fname)
    orb_beta = load_orbs_beta(fname)
    guess_core_hamiltonian(olp, kin + na, orb_alpha, orb_beta)
    orb_alpha.occupations[:3] = 1
    orb_beta.occupations[:] = 0
    dms = [orb_alpha.to_dm(), orb_beta.to_dm()]

    # converge scf and check the spins
    scf_solver = ODASCFSolver(1e-6)  # On some machines, 1e-8 does not work.
    scf_solver(ham, olp, occ_model, *dms)
    assert scf_solver.error(ham, olp, *dms) < scf_solver.threshold
    assert abs(np.einsum('ab,ba', olp, dms[0]) - 2) < 1e-10
    assert abs(np.einsum('ab,ba', olp, dms[1]) - 1) < 1e-10


def test_check_dm():
    # create random orthogonal vectors
    v1 = np.random.uniform(1, 2, 2)
    v1 /= np.linalg.norm(v1)
    v2 = np.array([-v1[1], v1[0]])
    v = np.array([v1, v2]).T

    olp = np.identity(2)

    op1 = np.dot(v * [-0.1, 0.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp)
    op1 = np.dot(v * [0.1, 1.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp)
    op1 = np.dot(v * [-0.1, 1.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp)
    op1 = np.dot(v * [0.1, 0.5], v.T)
    check_dm(op1, olp)
