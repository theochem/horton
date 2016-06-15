# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
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
from nose.plugins.attrib import attr

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.meanfield.test.common import check_hf_cs_hf, check_lih_os_hf, \
    check_water_cs_hfs, check_n2_cs_hfs, check_h3_os_hfs, check_h3_os_pbe, \
    check_co_cs_pbe, check_vanadium_sc_hf


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


def test_vanadium_sc_hf():
    with assert_raises(NoSCFConvergence):
        check_vanadium_sc_hf(ODASCFSolver(threshold=1e-10, maxiter=10))


def test_find_min_cubic():
    from horton.meanfield.scf_oda import find_min_cubic
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
    from horton.meanfield.scf_oda import find_min_quadratic
    assert find_min_quadratic(0.0, -0.7) == 1.0
    assert abs(find_min_quadratic(-3.0, 2.8) - 0.51724137931) < 1e-8
    assert abs(find_min_quadratic(-0.2, 0.1) - 0.666666666667) < 1e-8
    assert find_min_quadratic(-1.0, 1.0) == 0.5
    assert find_min_quadratic(-1.0, -1.0) == 1.0
    assert find_min_quadratic(1.0, 1.0) == 0.0
    assert find_min_quadratic(1.0, -0.5) == 0.0
    assert find_min_quadratic(1.0, -1.5) == 1.0


def test_aufbau_spin():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    mol = IOData.from_file(fn_fchk)
    occ_model = AufbauSpinOccModel(3)

    olp = mol.obasis.compute_overlap(mol.lf)
    kin = mol.obasis.compute_kinetic(mol.lf)
    na = mol.obasis.compute_nuclear_attraction(mol.coordinates, mol.pseudo_numbers, mol.lf)
    er = mol.obasis.compute_electron_repulsion(mol.lf)
    terms = [
        UTwoIndexTerm(kin, 'kin'),
        UDirectTerm(er, 'hartree'),
        UExchangeTerm(er,'x_hf'),
        UTwoIndexTerm(na, 'ne'),
    ]
    ham = UEffHam(terms)

    # Construct an initial state with instable spin polarization
    guess_core_hamiltonian(olp, kin, na, mol.exp_alpha, mol.exp_beta)
    mol.exp_alpha.occupations[:3] = 1
    mol.exp_beta.occupations[:] = 0
    dms = [mol.exp_alpha.to_dm(), mol.exp_beta.to_dm()]

    # converge scf and check the spins
    scf_solver = ODASCFSolver(1e-6) # On some machines, 1e-8 does not work.
    scf_solver(ham, mol.lf, olp, occ_model, *dms)
    assert scf_solver.error(ham, mol.lf, olp, *dms) < scf_solver.threshold
    assert abs(olp.contract_two('ab,ba', dms[0]) - 2) < 1e-10
    assert abs(olp.contract_two('ab,ba', dms[1]) - 1) < 1e-10



def test_check_dm():
    # create random orthogonal vectors
    v1 = np.random.uniform(1, 2, 2)
    v1 /= np.linalg.norm(v1)
    v2 = np.array([-v1[1], v1[0]])
    v = np.array([v1, v2]).T

    lf = DenseLinalgFactory(2)
    olp = lf.create_two_index()
    olp._array = np.identity(2)

    op1 = lf.create_two_index()
    op1._array = np.dot(v*[-0.1, 0.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp, lf)
    op1._array = np.dot(v*[0.1, 1.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp, lf)
    op1._array = np.dot(v*[-0.1, 1.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp, lf)
    op1._array = np.dot(v*[0.1, 0.5], v.T)
    check_dm(op1, olp, lf)
