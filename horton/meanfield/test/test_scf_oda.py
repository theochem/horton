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
#pylint: skip-file


import numpy as np
from nose.tools import assert_raises
from horton import *
from horton.meanfield.test.common import check_scf_hf_cs_hf, check_scf_water_cs_hfs


def test_scf_oda_cs_hf():
    check_scf_hf_cs_hf(SCFWrapper('oda', threshold=1e-7))


def test_scf_oda_cs_hfs():
    check_scf_water_cs_hfs(SCFWrapper('oda', threshold=1e-6))


def test_scf_oda_water_hf_321g():
    fn_fchk = context.get_fn('test/water_hfs_321g.fchk')
    sys = System.from_file(fn_fchk)
    scf_cache = Cache()
    ham = Hamiltonian(sys, scf_cache, [HartreeFockExchange(scf_cache, sys.lf, sys.wfn,
                                           sys.get_electron_repulsion())])

    # test continuation of interupted scf_oda
    guess_hamiltonian_core(sys)
    e0 = ham.compute()
    assert convergence_error_eigen(ham) > 1e-5
    with assert_raises(NoSCFConvergence):
        converge_scf_oda(ham, threshold=1e-2, maxiter=3)
    assert 'exp_alpha' in sys.wfn._cache
    e1 = ham.compute()
    with assert_raises(NoSCFConvergence):
        converge_scf_oda(ham, threshold=1e-2, maxiter=3)
    e2 = ham.compute()
    assert e1 < e0
    assert e2 < e1


def test_scf_oda_lih_hfs_321g():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    grid = BeckeMolGrid(sys.coordinates, sys.numbers, sys.pseudo_numbers, random_rotate=False)
    scf_cache = Cache()
    ham = Hamiltonian(sys, scf_cache, [Hartree(scf_cache, sys.lf, sys.wfn,
                                           sys.get_electron_repulsion()),
                            DiracExchange(scf_cache, sys.lf, sys.wfn,
                                           sys.get_electron_repulsion())],
                      grid)

    # test continuation of interupted scf_oda
    guess_hamiltonian_core(sys)
    e0 = ham.compute()
    assert convergence_error_eigen(ham) > 1e-5
    with assert_raises(NoSCFConvergence):
        converge_scf_oda(ham, threshold=1e-8, maxiter=3)
    assert 'exp_alpha' in sys.wfn._cache
    e1 = ham.compute()
    with assert_raises(NoSCFConvergence):
        converge_scf_oda(ham, threshold=1e-8, maxiter=3)
    e2 = ham.compute()
    assert e1 < e0
    assert e2 < e1

    # continue till convergence
    converge_scf_oda(ham, threshold=1e-8)


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


def test_scf_oda_aufbau_spin():
    fn_fchk = context.get_fn('test/li_h_3-21G_hf_g09.fchk')
    sys = System.from_file(fn_fchk)
    sys.wfn.occ_model = AufbauSpinOccModel(3)

    guess_hamiltonian_core(sys)
    scf_cache = Cache()
    ham = Hamiltonian(sys, scf_cache, [HartreeFockExchange(scf_cache, sys.lf, sys.wfn,
                                           sys.get_electron_repulsion())])
    converge_scf_oda(ham)


def test_check_dm():
    # create random orthogonal vectors
    v1 = np.random.uniform(1, 2, 2)
    v1 /= np.linalg.norm(v1)
    v2 = np.array([-v1[1], v1[0]])
    v = np.array([v1, v2]).T

    lf = DenseLinalgFactory(2)
    olp = lf.create_one_body()
    olp._array = np.identity(2)

    op1 = lf.create_one_body()
    op1._array = np.dot(v*[-0.1, 0.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp, lf, 'foo')
    op1._array = np.dot(v*[0.1, 1.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp, lf, 'foo')
    op1._array = np.dot(v*[-0.1, 1.5], v.T)
    with assert_raises(ValueError):
        check_dm(op1, olp, lf, 'foo')
    op1._array = np.dot(v*[0.1, 0.5], v.T)
    check_dm(op1, olp, lf, 'foo')
