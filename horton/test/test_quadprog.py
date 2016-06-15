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
from horton.quadprog import FeasibilityError, BoundedError, ConvergenceError, \
    _counter_to_free, find_1d_root, solve_safe, diagonal_form, \
    solve_constrained, solve_radius


def test_counter_to_free():
    free = np.zeros(5, bool)
    _counter_to_free(5, free)
    assert (free == [True, False, True, False, False]).all()
    _counter_to_free(15, free)
    assert (free == [True, True, True, True, False]).all()
    _counter_to_free(31, free)
    assert (free == [True, True, True, True, True]).all()
    _counter_to_free(16, free)
    assert (free == [False, False, False, False, True]).all()


def test_find_1d_root():
    # Only consider cases where the derivative of the function at the root is
    # at least 1.0. (See 'this is why'.)
    cases = [
        (-0.5, 0.1, np.sin, 0.0),
        (-10.0, 10.0, (lambda x: 2*(x-5)), 5.0),
        (-1.0, 2.0, (lambda x: np.exp(x)-1), 0.0),
        (-1.0, 2.0, (lambda x: np.exp(x)-2), np.log(2.0)),
        (0.0, 3.0, np.cos, np.pi/2),
    ]
    eps = 1e-5
    for x0, x2, fn, solution in cases:
        x1, y1 = find_1d_root(fn, (x0, fn(x0)), (x2, fn(x2)), eps)
        assert abs(y1) < eps
        assert abs(x1-solution) < eps # <--- this is why


def test_solve_safe():
    a = np.diag([1.0, 2.0, 1.0])
    b = np.array([2.0, 4.0, 1.0])
    assert (solve_safe(a, b) == [2.0, 2.0, 1.0]).all()
    a = np.diag([1.0, 2.0, 0.0])
    b = np.array([2.0, 4.0, 1.0])
    assert (solve_safe(a, b) == [2.0, 2.0, 0.0]).all()


def get_random_problem(nx=8, nl=None, posdef=False, epscn=1e-4):
    '''Generate a random constrainted quadratic problem

       **Optional arguments:**

       nx
            The number of unknowns

       nl
            The number of constraints. When not given, a random number from 0 to
            nx/2.

       posdef
            When set to True, a is positive definite.

       epscn
            The maximum condition number for a and r.
    '''
    if nl is None:
        nl = np.random.randint(0, nx/2)
    while True:
        # construct a
        a = np.random.normal(0, 1, (nx, nx))
        if posdef:
            a = np.dot(a, a.T)
        # resymmetrize in all cases because some blas implementations do not return
        # an exactly symmetric a times a^T.
        a = 0.5*(a + a.T)
        # avoid too ill-conditioned problems (rarely happens)
        absevals = abs(np.linalg.eigvalsh(a))
        if absevals.min() < epscn*absevals.max():
            continue
        b = np.random.normal(0, 1, nx)
        r = np.random.normal(0, 1, (nl, nx))
        if nl > 0:
            # avoid too ill-conditioned problems (rarely happens)
            abssingvals = abs(np.linalg.svd(r, full_matrices=False)[1])
            if abssingvals.min() < epscn*abssingvals.max():
                continue
        s = np.random.normal(0, 1, nl)
        return a, b, r, s


def test_diagonal_form_nl0():
    for i in xrange(100):
        # nl = 0
        a, b, r, s = get_random_problem(nl=0)
        x0, basis, evals, evecs, b_diag = diagonal_form(a, b, r, s)
        x1 = np.dot(evecs, b_diag/evals)
        x2 = solve_constrained(a, b, r, s)
        assert abs(x1 - x2).max() < 1e-8*abs(x1).max()
        assert x0 is None
        assert basis is None
        assert abs(np.dot(evecs*evals, evecs.T) - a).max() < 1e-10
        assert abs(np.dot(evecs, b_diag) - b).max() < 1e-10


def test_diagonal_form_nl4():
    for i in xrange(100):
        # nl = 5
        a, b, r, s = get_random_problem(nl=5)
        x0, basis, evals, evecs, b_diag = diagonal_form(a, b, r, s)
        assert abs(np.dot(r, x0) - s).max() < 1e-10
        assert abs(np.dot(r, basis.T)).max() < 1e-10
        assert evals.shape == (3,)
        assert evecs.shape == (3, 3)
        assert b_diag.shape == (3,)
        x3 = x0 + np.dot(basis.T, np.dot(evecs, b_diag/evals))
        x4 = solve_constrained(a, b, r, s)
        try:
            assert abs(x3 - x4).max() < 1e-8*abs(x3).max()
        except:
            print i, 5
            print np.linalg.eigvalsh(a)
            print np.linalg.svd(r)[1]
            print x3
            print x3 - x4
            print
            raise


def test_diagonal_form_nl8():
    for i in xrange(100):
        # nl = nx = 8
        a, b, r, s = get_random_problem(nl=8)
        x0, basis, evals, evecs, b_diag = diagonal_form(a, b, r, s)
        try:
            assert abs(np.dot(r, x0) - s).max() < 1e-10
        except:
            print i, 8, 'constraint'
            print np.linalg.svd(r)[1]
            print s
            print np.dot(r, x0) - s
            raise
        assert basis is None
        assert evals is None
        assert evecs is None
        assert b_diag is None
        x5 = solve_constrained(a, b, r, s)
        try:
            assert abs(x0 - x5).max() < 1e-8*abs(x0).max()
        except:
            print i, 8
            print np.linalg.eigvalsh(a)
            print np.linalg.svd(r)[1]
            print x0
            print x0 - x5
            print
            raise


def test_solve_radius_posdef_nl0():
    for i in xrange(100):
        a, b, r, s = get_random_problem(nl=0, posdef=True)
        radius = np.random.uniform(1e-5, 1e2)
        center = np.random.normal(0, 0.3*radius, len(b))
        x1 = solve_radius(a, b, center, radius, r, s)
        assert np.linalg.norm(x1 - center) < radius*(1 + 1e-3)
        x2 = solve_constrained(a, b, r, s)
        cost1 = np.dot(x1, 0.5*np.dot(a, x1) - b)
        cost2 = np.dot(x2, 0.5*np.dot(a, x2) - b)
        assert cost2 - cost1 < 1e-10
        if np.linalg.norm(x2 - center) < radius*(1 - 1e-3):
            assert abs(x1 - x2).max() < 1e-6*abs(x1).max()


def test_solve_radius_posdef_nl4():
    for i in xrange(100):
        a, b, r, s = get_random_problem(nl=4, posdef=True)
        radius = np.random.uniform(1e-5, 1e2)
        center = np.random.normal(0, 0.5*radius, len(b))
        try:
            x1 = solve_radius(a, b, center, radius, r, s)
        except FeasibilityError:
            continue
        assert np.linalg.norm(x1 - center) < radius*(1 + 1e-3)
        x2 = solve_constrained(a, b, r, s)
        cost1 = np.dot(x1, 0.5*np.dot(a, x1) - b)
        cost2 = np.dot(x2, 0.5*np.dot(a, x2) - b)
        assert cost2 - cost1 < 1e-10
        if np.linalg.norm(x2 - center) < radius*(1 - 1e-3):
            assert abs(x1 - x2).max() < 1e-6*abs(x1).max()


def test_solve_radius_posdef_nl8():
    for i in xrange(100):
        a, b, r, s = get_random_problem(nl=8, posdef=True)
        radius = np.random.uniform(1e-5, 1e2)
        center = np.random.normal(0, 0.5*radius, len(b))
        try:
            x1 = solve_radius(a, b, center, radius, r, s)
        except FeasibilityError:
            continue
        assert np.linalg.norm(x1 - center) < radius*(1 + 1e-3)
        x2 = solve_safe(r, s)
        assert abs(x1 - x2).max() < 1e-6*abs(x1).max()


def test_qps_infeasible():
    nx = 4
    nl = 2
    a = np.diag(np.random.uniform(2.0, 5.0, nx))
    b = np.random.normal(0.0, 1.0, nx)
    r = np.zeros((nl, nx))
    r[:,:2] = np.random.normal(0.0, 1.0, (nl, 2))
    s =  np.random.normal(0.0, 1.0, nl)
    qps = QPSolver(a, b, r, s)
    free = np.array([False, False, True, True])
    with assert_raises(FeasibilityError):
        qps.solve(free)


def test_qps_infeasible_radius():
    nx = 4
    nl = 2
    a = np.diag(np.random.uniform(2.0, 5.0, nx))
    b = np.random.normal(0.0, 1.0, nx)
    r = np.zeros((nl, nx))
    r[:,:2] = np.random.normal(0.0, 1.0, (nl, 2))
    s =  np.random.normal(0.0, 1.0, nl)
    center = np.random.normal(0, 1, nx)
    qps = QPSolver(a, b, r, s)
    free = np.array([False, False, True, True])
    with assert_raises(FeasibilityError):
        x = qps.solve_radius(free, center, 0.1)
        qps.check_feasible(x)


def get_random_feasible(qps, free=None, positive=False):
    if free is None:
        free = np.ones(qps.nx, bool)
    nfree = free.sum()
    for i in xrange(10):
        x_free = np.random.normal(0, 1, nfree)
        if qps.nl > 0:
            a_free, b_free, r_free, s_free = qps.get_free_problem(free)
            x_free += np.linalg.lstsq(r_free, s_free - np.dot(r_free, x_free))[0]
        if not positive or x_free.min() > 0:
            x = np.zeros(qps.nx)
            x[free] = x_free
            qps.check_feasible(x, free)
            return x


def get_random_free(qps):
    while True:
        free = np.random.randint(0, 2, qps.nx).astype(bool)
        if free.sum() > max(1, qps.nl):
            return free


def test_get_free_problem():
    for i0 in xrange(100):
        qps = QPSolver(*get_random_problem())
        for i1 in xrange(max(1, qps.nl), qps.nx+1):
            free = np.zeros(qps.nx, dtype=bool)
            free[:i1] = True
            a_free, b_free, r_free, s_free = qps.get_free_problem(free)
            assert a_free.shape == (i1, i1)
            assert b_free.shape == (i1,)
            if qps.nl > 0:
                assert r_free.shape == (qps.nl, i1)
                assert s_free.shape == (qps.nl,)
            else:
                assert r_free is None
                assert s_free is None


def test_expand_x():
    qps = QPSolver(*get_random_problem())
    for i0 in xrange(3, 9):
        free = np.zeros(qps.nx, dtype=bool)
        free[:i0] = True
        a_free, b_free, r_free, s_free = qps.get_free_problem(free)
        x_free = np.dot(a_free, b_free)
        x = qps._free_to_full(x_free)
        assert len(x) == qps.nx
        assert (x[i0:] == 0.0).all()
        assert (x[:i0] == x_free[:i0]).all()


def test_solve_qps():
    for i0 in xrange(100):
        qps = QPSolver(*get_random_problem())
        for i1 in xrange(10):
            if i1 == 0:
                free = np.ones(qps.nx, dtype=bool)
            else:
                free = get_random_free(qps)
            x = qps.solve(free)
            qps.check_feasible(x, free)


def test_solve_qps_radius_posdef():
    for i0 in xrange(100):
        qps = QPSolver(*get_random_problem(posdef=True))
        for i1 in xrange(10):
            free = get_random_free(qps)
            center = get_random_feasible(qps, free)
            if center is None:
                continue
            radius = np.random.uniform(0.1, 0.5)
            cost_center = qps.compute_cost(center)
            x = qps.solve_radius(free, center, radius)
            qps.check_feasible(x, free)
            cost_x = qps.compute_cost(x)
            assert cost_x < cost_center
            assert np.linalg.norm(x - center) <= radius*1.001


def test_solve_qps_radius():
    for i0 in xrange(100):
        qps = QPSolver(*get_random_problem())
        for i1 in xrange(10):
            free = get_random_free(qps)
            center = get_random_feasible(qps, free)
            if center is None:
                continue
            radius = np.random.uniform(0.1, 0.5)
            cost_center = qps.compute_cost(center)
            x = qps.solve_radius(free, center, radius)
            qps.check_feasible(x, free)
            cost_x = qps.compute_cost(x)
            assert cost_x < cost_center
            assert np.linalg.norm(x - center) <= radius*1.001


def test_brute_simple():
    a = np.identity(5)
    b = -2*np.ones(5)
    r = np.ones((1, 5))
    s = np.array([5])
    qps = QPSolver(a, b, r, s)
    cost, x = qps.find_brute()
    assert abs(x - np.ones(5)).max() < 1e-10
    qps.check_solution(x)


@attr('slow')
def test_brute_local_posdef():
    for counter in xrange(100):
        # A large eps is used because some random problems are very ill-behaved.
        qps = QPSolver(*get_random_problem(nx=6, posdef=True), eps=1e-6)

        try:
            cost0, x0 = qps.find_brute()
            feasible = True
        except FeasibilityError:
            feasible = False

        if feasible:
            try:
                qps.check_solution(x0)
            except:
                print 'problem with brute'
                qps.log()
                raise

            try:
                cost1, x1 = qps.find_local(x0, 100.0)
                qps.check_solution(x1)
                assert abs(x0 - x1).max() < qps.eps
            except:
                print 'problem with local from solution'
                qps.log(x0)
                raise

            xr = get_random_feasible(qps, positive=True)
            if xr is None:
                continue
            try:
                cost2, x2 = qps.find_local(xr, 100.0)
                qps.check_solution(x2)
                assert abs(x0 - x2).max() < qps.eps
            except:
                print 'problem with local from random'
                qps.log(xr)
                raise


@attr('slow')
def test_brute_local():
    for counter in xrange(100):
        # A large eps is used because some random problems are very ill-behaved.
        qps = QPSolver(*get_random_problem(nx=6), eps=1e-6)
        qps.log()

        try:
            cost0, x0 = qps.find_brute()
            feasible = True
        except (FeasibilityError, BoundedError), e:
            feasible = False

        if feasible:
            try:
                qps.check_solution(x0)
            except:
                print 'problem with brute'
                qps.log()
                raise

            try:
                cost1, x1 = qps.find_local(x0, 100.0)
                qps.check_solution(x1)
                assert abs(x0 - x1).max() < qps.eps
            except:
                print 'problem with local from solution'
                qps.log(x0)
                raise

            xr = get_random_feasible(qps, positive=True)
            qps.log(xr)
            if xr is None:
                continue
            try:
                cost2, x2 = qps.find_local(xr, 100.0)
                qps.check_solution(x2)
                assert cost0 - cost2 < qps.eps
            except ConvergenceError:
                # this may happen with unbounded problems
                continue
            except:
                print 'problem with local from random'
                qps.log(xr)
                raise


def test_brute_case1():
    qps = QPSolver(
        np.array([
            [5.8818949999422374, -0.77419795370541555, 0.39341617187739442, 2.0623472342526616, 0.70305696723217903, 0.43927082926019057],
            [-0.77419795370541555, 1.2417514864801085, -0.87970058018596486, -2.2750064483132406, -0.53127538666014551, 0.54723718046867531],
            [0.39341617187739442, -0.87970058018596486, 8.1424259636396421, 2.0743354833933023, -0.8451551155934458, -0.87608149252282963],
            [2.0623472342526616, -2.2750064483132406, 2.0743354833933023, 6.0910659702544807, -1.0972631644140005, -0.23028504167149105],
            [0.70305696723217903, -0.53127538666014551, -0.8451551155934458, -1.0972631644140005, 5.565314789423871, -3.7575632324834198],
            [0.43927082926019057, 0.54723718046867531, -0.87608149252282963, -0.23028504167149105, -3.7575632324834198, 5.2499593294744216],
        ]),
        np.array([0.79794082940662259, 1.7133727997077115, 0.52660856271888534, 0.87508447140496493, 1.4399900959090817, 1.890635603184891]),
        np.array([
            [-0.3066400060217162, 0.4347344395513979, -2.3801164580913214, -0.090650765867039226, 0.086345721171055323, 1.1449625586344025],
            [-0.279502458226559, 0.59702302554570619, 0.21989715465365572, -0.11013881571319613, 0.16478641607527697, 0.93466512224244058],
        ]),
        np.array([-0.99746332627189982, 1.1089528830442121]),
    )
    cost, x = qps.find_brute()
    qps.check_solution(x)


def test_brute_case2():
    qps = QPSolver(
        np.array([
            [7.7777678491839488, -4.1394215100441691, -0.93496478768154112, -2.5636063412433505, 0.4152351107141245, 2.1325362349317176],
            [-4.1394215100441691, 5.3966791477501044, -0.078655286246388428, 0.84760054283629427, 0.85251884466077499, -4.1538822573446677],
            [-0.93496478768154112, -0.078655286246388428, 3.2776150136403599, 0.66562278015287957, -2.01031234900845, -1.3170956617298368],
            [-2.5636063412433505, 0.84760054283629427, 0.66562278015287957, 8.0220111474488096, -2.8503515897516283, -1.2294009495541263],
            [0.4152351107141245, 0.85251884466077499, -2.01031234900845, -2.8503515897516283, 2.3078909975036965, 0.51158004287230041],
            [2.1325362349317176, -4.1538822573446677, -1.3170956617298368, -1.2294009495541263, 0.51158004287230041, 4.753207067427045],
        ]),
        np.array([0.084122796400789776, 0.13248005552515144, -0.011337219032290116, -0.32808318842735468, -0.1942335218652631, 0.42946358838556126]),
        np.array([
            [0.058947500030390558, 2.8343622100824701, 1.3178663807154574, 0.63164184657317068, 0.83782836334183131, -0.051926442369153891],
            [-0.31073645606613759, -0.50172261649450467, 0.1826071399900123, 1.5181979800166607, -0.67838760926895592, -0.40545902742915829],
        ]),
        np.array([-0.57685125471751075, -1.2399862660238066]),
    )
    with assert_raises(FeasibilityError):
        cost, x = qps.find_brute()


def test_brute_case3():
    qps = QPSolver(
        np.array([
            [10.949643502677578, -6.4848706832640168, -2.2855976904588635, 6.2603152305513392, -1.0445724421867459, 2.9962572283103555],
            [-6.4848706832640168, 8.3170212543577193, 0.27469715271157841, -4.6215940087176399, 2.712716602561124, -4.8477069972400084],
            [-2.2855976904588635, 0.27469715271157841, 14.769272405279011, -0.051505059995286595, -0.49106152710130235, -2.1282377916520616],
            [6.2603152305513392, -4.6215940087176399, -0.051505059995286595, 6.772933209165914, 0.31613657279823015, 4.2766554080787502],
            [-1.0445724421867459, 2.712716602561124, -0.49106152710130235, 0.31613657279823015, 3.741204058931006, -1.7853429804383523],
            [2.9962572283103555, -4.8477069972400084, -2.1282377916520616, 4.2766554080787502, -1.7853429804383523, 6.9502814568206386],
        ]),
        np.array([0.81811916699197473, -0.2839108577501086, -0.52580285068573052, 1.1168945191958584, 0.51473677332239764, 0.40132813088081393]),
        np.array([
            [0.60693853800344977, -1.602707177915851, -0.045461208556361564, -1.3574337407982444, -1.0194532203245703, -0.11007054336981964],
            [0.41595858414246062, 0.20880157188047851, -1.1189280276138425, -1.3021514511281107, 0.10768134383771397, 1.7556828215121263],
        ]),
        np.array([0.51812694007155591, -0.82287154751311653]),
    )
    cost, x = qps.find_brute()
    qps.check_solution(x)


def test_brute_case4():
    qps = QPSolver(
        np.array([
            [0.7990546151426382, 0.60159970689434261, 0.25038104422335095, 0.9107186394321235, -0.057596715716904341, 0.36415173883112567],
            [0.60159970689434261, 0.071406318978636474, -0.18682517604341942, 0.1408876305178936, -1.4706046610777006, 0.95698263756937496],
            [0.25038104422335095, -0.18682517604341942, -0.76023253647318112, -0.65971813975561078, 0.091128703324063975, 0.72341372601515963],
            [0.9107186394321235, 0.1408876305178936, -0.65971813975561078, -1.3564464302489747, 0.14486671746488211, 0.37478117023135221],
            [-0.057596715716904341, -1.4706046610777006, 0.091128703324063975, 0.14486671746488211, -0.44189735537298425, 0.77868082771555558],
            [0.36415173883112567, 0.95698263756937496, 0.72341372601515963, 0.37478117023135221, 0.77868082771555558, 0.15083792086070447],
        ]),
        np.array([-2.1492020469280524, 1.5260662967668095, 0.91429181277522631, -0.26822566090243904, -0.34729383842661415, 0.67908749141704505]),
        np.array([
            [-0.41015693790806101, -0.3029115511740868, -0.68741014841460135, 0.16295599958248516, -0.8286427429967762, -0.22428931788921946],
        ]),
        np.array([-0.75113559150400466]),
    )
    with assert_raises(BoundedError):
        cost, x = qps.find_brute()


def test_brute_case5():
    qps = QPSolver(
        np.array([
            [0, 0.18406347881903429],
            [0.18406347881903429, 0],
        ]),
        np.array([97.843415950141335, 98.560421251173764]),
        np.array([
            [1.0, 1.0],
        ]),
        np.array([1.0]),
    )
    assert qps.compute_cost(np.array([1.0, 0.0])) > qps.compute_cost(np.array([0.0, 1.0]))
    cost, x = qps.find_brute()
    qps.check_solution(x)


def test_brute_case6():
    # all fixed is not a good idea ...
    qps = QPSolver(
        np.array([
            [-1.6170622392405418, -1.1110660450798633, 0.15534256626761522, -0.035521902704618108, -0.7461985258486824, 0.10251918427316184],
            [-1.1110660450798633, -0.34288398900411249, -0.87112101192198466, -1.0723409286884387, 0.77890683895688473, 0.64749152585353753],
            [0.15534256626761522, -0.87112101192198466, -0.47456930397068336, 1.0664188304998987, 0.52494441816706661, -1.1798151543302824],
            [-0.035521902704618108, -1.0723409286884387, 1.0664188304998987, -0.062964411407941484, -0.2032160752170718, 1.7280732016323954],
            [-0.7461985258486824, 0.77890683895688473, 0.52494441816706661, -0.2032160752170718, -1.6672524270549614, -0.0079015689549524204],
            [0.10251918427316184, 0.64749152585353753, -1.1798151543302824, 1.7280732016323954, -0.0079015689549524204, -2.4093790812393414],
        ]),
        np.array([-1.628027986608761, -1.680081567157873, 0.12815213529315533, -0.38201939028605525, 0.77151603852902884, 0.50875081807696942]),
        np.zeros((0,6)),
        np.array([]),
    )
    with assert_raises(BoundedError):
        cost, x = qps.find_brute()


def test_brute_case7():
    qps = QPSolver(
        np.array([
            [-0, -2.3234747459355276e-12, -2.8634872251132037e-12, -2.2701840407535201e-12],
            [-2.3234747459355276e-12, -0, -3.907985046680551e-14, -0],
            [-2.8634872251132037e-12, -3.907985046680551e-14, -0, -3.1974423109204508e-14],
            [-2.2701840407535201e-12, -0, -3.1974423109204508e-14, -0],
        ]),
        np.array([97.540763972605916, 97.540763972610776, 97.540763972610804, 97.540763972610776]),
        np.array([
            [1.0, 1.0, 1.0, 1.0],
        ]),
        np.array([1.0]),
    )
    cost, x = qps.find_brute()
    qps.check_solution(x)


def test_brute_case8():
    # find_brute returns a solution that is not exactly on the equality constraints.
    qps = QPSolver(
        np.array([
            [-0.70841731236733996, -1.2971293403402491, 0.90762745129580269, 0.086526845993940282, 0.83018921812606727, 0.31381838847861343],
            [-1.2971293403402491, -1.9408646843381616, -0.89540919308959155, -0.024479335288883169, 1.0775985944280677, -0.52650657733203432],
            [0.90762745129580269, -0.89540919308959155, 0.40612490363321968, 0.014668804820934134, 0.5370600831832697, -0.39964394635396933],
            [0.086526845993940282, -0.024479335288883169, 0.014668804820934134, 0.015929091021452374, -0.79212879441047745, 0.016890690484653692],
            [0.83018921812606727, 1.0775985944280677, 0.5370600831832697, -0.79212879441047745, -0.88786447750589015, -0.034978317332826636],
            [0.31381838847861343, -0.52650657733203432, -0.39964394635396933, 0.016890690484653692, -0.034978317332826636, -0.18343755938048792],
        ]),
        np.array([-0.97667398736166677, -0.26674955274778495, 1.3971837428664151, -0.7548208542900311, 0.20143172017397515, 0.70870382939012788]),
        np.array([
            [0.23201895073676565, 1.2815806707407724, -0.48070981450845679, -0.086198330259237302, -0.24520586286913254, -0.17363911537396109],
            [-0.68428995101920931, -1.0316773608304244, -0.17545371137719318, -0.2729594433963799, 0.19735403185101402, -0.97249697065508856],
        ]),
        np.array([0.41079612937870091, -1.6343480803544916]),
        eps=1e-6,
    )
    cost, x = qps.find_brute()
    qps.check_solution(x)


def test_brute_case9():
    # mismatch equality constraints (posdef)
    qps = QPSolver(
        np.array([
            [4.3206954954195442, 1.0995193491529669, -1.8495021618983287, 0.72617849172585236, -1.7933578310574625, -3.3858794413657338],
            [1.0995193491529669, 7.4747802231556353, -7.4921796351274246, -1.5044955980863084, -2.9385575178317112, -2.5930622344533512],
            [-1.8495021618983287, -7.4921796351274246, 11.271691407429426, -0.13747973098167873, 3.0801710797885811, 2.9305864693055015],
            [0.72617849172585236, -1.5044955980863084, -0.13747973098167873, 7.1223220952346562, 6.3309090398866426, 0.89431284768332786],
            [-1.7933578310574625, -2.9385575178317112, 3.0801710797885811, 6.3309090398866426, 10.48601433755171, 2.3729785411834676],
            [-3.3858794413657338, -2.5930622344533512, 2.9305864693055015, 0.89431284768332786, 2.3729785411834676, 3.6504940434519155],
        ]),
        np.array([-0.20354699990469935, 1.8853961326540016, -0.7319189254058831, 0.065576550436379791, -0.94479064660325729, 1.0746348777770884]),
        np.array([
            [1.2485461200060213, 1.2175464308211492, 1.0092841327024706, 2.7899706809845943, -0.92878294027396968, -0.35937248238222808],
            [-0.31042069923722482, 1.3203168967639625, 0.71261924983670422, 1.199582577635893, 0.2304738785590868, 0.15518540117148097],
        ]),
        np.array([0.40260474485255499, -0.48930086847135434]),
        eps=1e-6,
    )
    cost, x = qps.find_brute()
    qps.check_solution(x)


def test_local_case1():
    qps = QPSolver(
        np.array([
            [12.98721885864512, -4.6588586635198244, -0.91546482655448469, -5.9789922145586925, 2.7503595102379421, 1.5317082550812957],
            [-4.6588586635198244, 7.0326437030982909, -1.2068162572718182, 4.1152385083361009, -0.486080361964501, -0.92288125639729679],
            [-0.91546482655448469, -1.2068162572718182, 7.0824548119378132, -3.4096249540917754, -0.9217331851519901, -1.5005654709949769],
            [-5.9789922145586925, 4.1152385083361009, -3.4096249540917754, 11.89628880624144, -3.6316543061999216, 0.42944160122592212],
            [2.7503595102379421, -0.486080361964501, -0.9217331851519901, -3.6316543061999216, 3.3531338144290426, 0.49086672463256498],
            [1.5317082550812957, -0.92288125639729679, -1.5005654709949769, 0.42944160122592212, 0.49086672463256498, 0.86607069729952701],
        ]),
        np.array([-1.7640856679361303, -0.57087943458650181, -0.13374131333031397, 0.34522781648137341, 0.068577833332047478, -0.11172137814018418]),
        np.array([
            [0.43101465265031758, 0.086243316821451532, 0.056467524174045318, -0.83513643817703564, -0.95396057797839307, -1.1277118607292189],
            [-0.063668029637385076, 1.2248051747709392, 0.17965412995933414, -0.55365467576520333, -0.42904406085688018, -1.9683318875692319],
        ]),
        np.array([0.65794380142477105, 0.038165380592388352]),
    )
    guess = np.array([1.4709241946911222, 0.058506557259324535, 0.33484947848787516, 0.0, 0.0, 0.0])
    cost, x = qps.find_local(guess, 1.0)
    qps.check_solution(x)


def test_local_case2():
    # problem with new_cost < cost
    qps = QPSolver(
        np.array([
            [7.4629736594322331, -0.78748447995959669, 0.39350017508027152, -0.68307208660912533, 4.9135859685012173, 1.3191440233541152],
            [-0.78748447995959669, 2.8455612103600143, -2.3777176627456003, -2.0647887688885893, -0.8494420668115783, 1.0362029379461384],
            [0.39350017508027152, -2.3777176627456003, 7.1282320417634413, 1.1073706764518034, 1.6754999335349345, 0.24762190846137289],
            [-0.68307208660912533, -2.0647887688885893, 1.1073706764518034, 2.4252041096631261, -0.69785432183912288, -0.55618812054465072],
            [4.9135859685012173, -0.8494420668115783, 1.6754999335349345, -0.69785432183912288, 4.4961443889941322, 0.56183256135502446],
            [1.3191440233541152, 1.0362029379461384, 0.24762190846137289, -0.55618812054465072, 0.56183256135502446, 1.9755149613389364],
        ]),
        np.array([-0.036326161624374083, 0.93442606896402636, 2.3465206569805046, 0.60013364713670336, -0.065360649313220678, 0.74605936179230847]),
        np.zeros((0,6)),
        np.array([]),
    )
    guess = np.array([0.060781020807135232, 0.094832047408032777, 1.1736098312337082, 0.27679006326937744, 0.08915412896212016, 2.7251320294146675])
    cost, x = qps.find_local(guess, 1.0)
    qps.check_solution(x)


def test_local_case3():
    # problem with tmin > 0
    qps = QPSolver(
        np.array([
            [3.0416282731277371, -1.4216164332936789, 0.46036692207941887, 0.8646633606722256, 0.37224511445850084, 0.49853885743325299],
            [-1.4216164332936789, 10.48010471548946, -5.0724948663530656, 0.51685701382721216, -5.7757019664518729, -2.654913315245504],
            [0.46036692207941887, -5.0724948663530656, 4.6454659562180787, -2.8459087108315506, 2.2400328467875963, 3.5852631057103288],
            [0.8646633606722256, 0.51685701382721216, -2.8459087108315506, 6.9625681947633717, 0.0065222059726193216, -4.072767567818075],
            [0.37224511445850084, -5.7757019664518729, 2.2400328467875963, 0.0065222059726193216, 7.0132899829952953, -0.99707593216694956],
            [0.49853885743325299, -2.654913315245504, 3.5852631057103288, -4.072767567818075, -0.99707593216694956, 4.6945205762054307],
        ]),
        np.array([-0.45709100400275965, -1.8550340629069841, -0.22415042613116959, 0.42530295591391204, 0.99077132382609068, 0.44684458290517975]),
        np.array([
            [-0.32925296164525752, 0.34285843523483472, -1.1876007007262734, -0.34405996689960999, -0.48748723186183124, -0.84790951916486579],
            [0.34748864320332279, -1.4306241151049646, 0.48663092329413032, 0.36880547559765398, 0.5830781797631559, 0.56629405720887827],
        ]),
        np.array([-0.8042070619771764, -0.8113771733846965]),
    )
    guess = np.array([0.33919202323040087, 1.4357688663259185, 0.128273057926652, 0.84509497913350029, 0.99165656530820734, 0.30459867095075666])
    cost, x = qps.find_local(guess, 1.0)
    qps.check_solution(x)


def test_local_case4():
    # Impossible to converge
    qps = QPSolver(
        np.array([
            [-0.41833209771758612, 0.062887697148893529, -0.016088048200215366, -0.050445141412105376, 0.28111333503256047, 0.34422748081701854],
            [0.062887697148893529, -0.37785375633297413, 0.95384534852986791, -0.17966631075461084, -0.54690073408327788, 0.74422105092988078],
            [-0.016088048200215366, 0.95384534852986791, -0.93407522326505077, 0.46787002711838876, 0.87574341703436942, -0.74308895418662135],
            [-0.050445141412105376, -0.17966631075461084, 0.46787002711838876, 0.095602896288350628, -0.43330607692880641, -0.09673966777597498],
            [0.28111333503256047, -0.54690073408327788, 0.87574341703436942, -0.43330607692880641, 0.08109348599752926, 1.8753235082158164],
            [0.34422748081701854, 0.74422105092988078, -0.74308895418662135, -0.09673966777597498, 1.8753235082158164, -0.29575516753680142],
        ]),
        np.array([-1.369335961961399, 1.167500700090857, 0.21414477811997978, 2.5649100860678695, -1.9806821349849604, -1.5754238007215169]),
        np.array([
            [-0.41956342254107903, 0.71674035950201276, -1.5770905037800305, 0.43181073722954405, -0.97223306449577329, -0.62828793519627313],
        ]),
        np.array([-0.33651216328153544]),
    )
    guess = np.array([1.4532630642344062, 0.72265262030217081, 0.065281886690390367, 0.1839421913168558, 0.20535723808514028, 0.03429655200129661])
    with assert_raises(ConvergenceError):
        qps.find_local(guess, 1.0)


def test_local_case5():
    # slow convergence
    qps = QPSolver(
        np.array([
            [0.79079907934182903, 0.011808066578099186, -1.0980998280784657, -0.68897546615676319, -0.10555819609322133, -0.0094175398403186783],
            [0.011808066578099186, -1.2570947595949307, -0.25176742901500682, 0.43452314945243631, -0.51301954393280846, -0.31209650217715923],
            [-1.0980998280784657, -0.25176742901500682, 1.6781911482939742, -0.70553038981825766, 0.54040995405217274, 0.10204140147104501],
            [-0.68897546615676319, 0.43452314945243631, -0.70553038981825766, 0.12542537761324651, 0.40866389380559537, 0.18486598005512384],
            [-0.10555819609322133, -0.51301954393280846, 0.54040995405217274, 0.40866389380559537, -0.92415498156459686, 0.95932603388371385],
            [-0.0094175398403186783, -0.31209650217715923, 0.10204140147104501, 0.18486598005512384, 0.95932603388371385, 0.23569100050886349],
        ]),
        np.array([-0.46333818345367495, 0.21257179073800472, 0.62688340581190316, 2.2623297911779567, 0.50819379925975028, -0.6943526670820328]),
        np.array([
            [-1.2657558954797123, -0.58071957697534726, 1.7286979572389156, 0.79253342400981885, -0.92638488825515308, 0.42858228723662439],
            [-0.91371314768106349, -0.33354701632644296, 0.29090276846263258, 0.025115507692647247, -1.4686929816358802, 0.12857290086890488],
        ]),
        np.array([2.7532334775271505, -0.48328722498961846]),
    )
    guess = np.array([0.17076836616570085, 0.11727673011988589, 1.1886537627251008, 1.675822839903278, 0.48066423537800507, 0.23286211070422069])
    cost, x = qps.find_local(guess, 1.0)
    qps.check_solution(x)


def test_local_case6():
    # positive definite case that has convergence problems.
    # the trust radius should not be too small... solved.
    qps = QPSolver(
        np.array([
            [3.049138030515488, -2.267339701331577, 1.6452178953594849, -1.1114025762214919, -1.4683391262333996, -2.8949938905463575],
            [-2.267339701331577, 4.0639013616185977, -0.85512669477123537, 0.16765289716466125, 0.54721401159266114, 2.3075404741986758],
            [1.6452178953594849, -0.85512669477123537, 4.4883269223774871, -5.4322963136631426, -2.5999282829704065, -4.2203649379170152],
            [-1.1114025762214919, 0.16765289716466125, -5.4322963136631426, 11.291537524605506, 2.1266925840689384, -0.15939236254485001],
            [-1.4683391262333996, 0.54721401159266114, -2.5999282829704065, 2.1266925840689384, 3.3080344785127451, 2.7570212908688965],
            [-2.8949938905463575, 2.3075404741986758, -4.2203649379170152, -0.15939236254485001, 2.7570212908688965, 10.413810265914657],
        ]),
        np.array([2.1051178688693808, 0.77034513234528068, -0.066413262422188549, -1.161467599974259, -1.145284458324215, 1.9115896887435808]),
        np.zeros((0,6)),
        np.array([]),
    )
    guess = np.array([1.1218006104444309, 0.055845280862766397, 0.5678162968457573, 0.64202055422478232, 0.79817332647621919, 0.4492255175376717])
    cost, x = qps.find_local(guess, 10.0)
    qps.check_solution(x)


def test_local_case7():
    # local search seems to beat global one?
    qps = QPSolver(
        np.array([
            [-0.70579902180480114, -0.4327726229988561, 0.24228757421240607, 0.51418623618564996, 1.0632780675093283, -1.2852213384501507],
            [-0.4327726229988561, -1.0084276903656637, 0.27249303018293985, 1.3760173009608194, 0.67570216250059412, 0.52160016838454792],
            [0.24228757421240607, 0.27249303018293985, -0.125476715180112, -0.21437197456564006, 0.92014307041811716, 0.77230951003769155],
            [0.51418623618564996, 1.3760173009608194, -0.21437197456564006, 0.11169101328138263, -0.58757562610076719, 0.54870518654715494],
            [1.0632780675093283, 0.67570216250059412, 0.92014307041811716, -0.58757562610076719, 1.3791708665173121, -0.44585251609595422],
            [-1.2852213384501507, 0.52160016838454792, 0.77230951003769155, 0.54870518654715494, -0.44585251609595422, -1.5140735372378349],
        ]),
        np.array([-2.2854963174099305, -0.70412551963528325, -0.20443019124351977, -1.750838749153953, 0.92664898509539495, -1.3205815413488067]),
        np.array([
            [-1.4834034545177326, 0.61696635160928115, 1.2902252770986995, -0.96826004067118876, -1.865873147282245, 2.1304980879954329],
        ]),
        np.array([-1.8919484511303895]),
    )
    cost0, x0 = qps.find_brute()
    qps.check_solution(x0)
    guess = np.array([0.66874113007080915, 0.38269689497818871, 0.43551620861147256, 1.7616957362635119, 0.87395834895248681, 0.76907570723573859])
    cost1, x1 = qps.find_local(guess, 10.0)
    qps.check_solution(x1)
    assert cost0 - cost1 < qps.eps


def test_local_case8():
    # This comes from an old implementation of EDIIS. The initial guess was
    # not constructed properly. Anyway, with proper initial ridge parameters,
    # this also works.
    qps = QPSolver(
        np.array([
            [-0, -0.80335249114943963, -0.4810486491314947, -0.089239349136249757, -0.15436494944003698, -0.35230777891029419, -0.28555315037639062, -0.19177912108670014, -0.21585294726538962, -0.26476652332778627],
            [-0.80335249114943963, -0, -0.040921462946748477, -0.38027796950886028, -0.2660150289073222, -0.098718781279508505, -0.13982199557364083, -0.22452799188664585, -0.19912180157206549, -0.15608175980106864],
            [-0.4810486491314947, -0.040921462946748477, -0, -0.17321382433037158, -0.099670821942954291, -0.014896131374744215, -0.031613318232015075, -0.075875344294711766, -0.061676589860272912, -0.039439435052052829],
            [-0.089239349136249757, -0.38027796950886028, -0.17321382433037158, -0, -0.010168990909519948, -0.092711943454681034, -0.059610136408235093, -0.02078877416074576, -0.029532988876905364, -0.049855163352130916],
            [-0.15436494944003698, -0.2660150289073222, -0.099670821942954291, -0.010168990909519948, -0, -0.041590858093009331, -0.020646513109756626, -0.0020336615488893983, -0.0051892404310116547, -0.015133796800792965],
            [-0.35230777891029419, -0.098718781279508505, -0.014896131374744215, -0.092711943454681034, -0.041590858093009331, -0, -0.0036392226383412662, -0.025764057353654835, -0.017631978012037308, -0.0066142161822391188],
            [-0.28555315037639062, -0.13982199557364083, -0.031613318232015075, -0.059610136408235093, -0.020646513109756626, -0.0036392226383412662, -0, -0.010041713436219624, -0.0052548807969650113, -0.00044679245518608468],
            [-0.19177912108670014, -0.22452799188664585, -0.075875344294711766, -0.02078877416074576, -0.0020336615488893983, -0.025764057353654835, -0.010041713436219624, -0, -0.00076870247709948103, -0.006272367939274659],
            [-0.21585294726538962, -0.19912180157206549, -0.061676589860272912, -0.029532988876905364, -0.0051892404310116547, -0.017631978012037308, -0.0052548807969650113, -0.00076870247709948103, -0, -0.0026493882418954229],
            [-0.26476652332778627, -0.15608175980106864, -0.039439435052052829, -0.049855163352130916, -0.015133796800792965, -0.0066142161822391188, -0.00044679245518608468, -0.006272367939274659, -0.0026493882418954229, -0],
        ]),
        np.array([96.74931669030984, 97.210304973486174, 97.2575693294621, 97.462260931861223, 97.472700776115133, 97.518784700439411, 97.5224733817101, 97.535582111324103, 97.536350310192589, 97.539372850065561]),
        np.array([
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ]),
        np.array([1.0]),
    )
    guess = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    cost, x = qps.find_local(guess, 1.0)
    qps.check_solution(x)


def test_local_case9():
    # yet another hard case for convergence
    qps = QPSolver(
        np.array([
            [-0, -241.45520097844334, -81.025811808215622, -66.646386842604926, -90.984366008742555, -52.004805821364386, -70.425468607943742, -57.366620069719914, -55.1877804833006, -53.705769566242907, -56.391096977161254],
            [-241.45520097844334, -0, -162.135903992414, -90.262677632643232, -70.276892079044018, -92.875705599218321, -72.28463739472835, -79.197990735263048, -82.586075985622472, -83.994287707438332, -78.796070159949906],
            [-81.025811808215622, -162.135903992414, -0, -87.28951414271836, -48.242358569534147, -28.742953020160229, -49.101136904866735, -32.105000819299548, -40.850665783027893, -30.057423682337586, -33.735044181468567],
            [-66.646386842604926, -90.262677632643232, -87.28951414271836, -0, -57.659521455279005, -24.69684419657078, -20.307639293415008, -21.441458536207293, -11.199312043260591, -19.974590362450328, -17.916540652877458],
            [-90.984366008742555, -70.276892079044018, -48.242358569534147, -57.659521455279005, -0, -27.948974737346248, -13.712164271081321, -14.218021097302355, -30.862156391808355, -19.790364207820602, -17.217970999822271],
            [-52.004805821364386, -92.875705599218321, -28.742953020160229, -24.69684419657078, -27.948974737346248, -0, -18.060952321337993, -4.025388940084369, -3.8482396016948144, -1.8363489948160066, -3.644582361475031],
            [-70.425468607943742, -72.28463739472835, -49.101136904866735, -20.307639293415008, -13.712164271081321, -18.060952321337993, -0, -5.93835013335962, -10.593591943371536, -8.6487068864670533, -5.644151180036971],
            [-57.366620069719914, -79.197990735263048, -32.105000819299548, -21.441458536207293, -14.218021097302355, -4.025388940084369, -5.93835013335962, -0, -4.241229635355424, -0.7847837230720387, -0.29680249521985047],
            [-55.1877804833006, -82.586075985622472, -40.850665783027893, -11.199312043260591, -30.862156391808355, -3.8482396016948144, -10.593591943371536, -4.241229635355424, -0, -2.1191260827439748, -2.4949817756629358],
            [-53.705769566242907, -83.994287707438332, -30.057423682337586, -19.974590362450328, -19.790364207820602, -1.8363489948160066, -8.6487068864670533, -0.7847837230720387, -2.1191260827439748, -0, -0.46855762948329982],
            [-56.391096977161254, -78.796070159949906, -33.735044181468567, -17.916540652877458, -17.217970999822271, -3.644582361475031, -5.644151180036971, -0.29680249521985047, -2.4949817756629358, -0.46855762948329982, -0],
        ]),
        np.array([219.88503892388073, 206.47122847047703, 264.95062224390711, 282.07819206858949, 285.59235059522575, 303.26503791957691, 298.61126236074369, 307.14088474324785, 304.05986518939471, 307.14000745089061, 307.31174519920341]),
        np.array([
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        ]),
        np.array([1.0]),
    )
    guess = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
    cost, x = qps.find_local(guess, 10.0)
    qps.check_solution(x)


def test_local_case10():
    # failed to converge (posdef)
    # This turned out to be a trust radius that was too small.
    qps = QPSolver(
        np.array([
            [3.6420954087755524, -3.0520276574272556, 0.8600212081041877, -0.27002779932504622, -1.4149411955120152, -1.7656114953692559],
            [-3.0520276574272556, 5.2980117879019053, -0.51863311310056681, -3.1078082667477691, 1.2965118123600228, 2.5974001710229837],
            [0.8600212081041877, -0.51863311310056681, 2.2859933272503223, 0.38632978005515767, 0.24074547614702529, -2.5864818874051134],
            [-0.27002779932504622, -3.1078082667477691, 0.38632978005515767, 4.5718542781367075, 0.31035781257497425, -2.2997148639134855],
            [-1.4149411955120152, 1.2965118123600228, 0.24074547614702529, 0.31035781257497425, 3.592559354700823, -0.32310839292740823],
            [-1.7656114953692559, 2.5974001710229837, -2.5864818874051134, -2.2997148639134855, -0.32310839292740823, 4.1635594436043855],
        ]),
        np.array([0.6827237762636329, 0.61262449970997035, -0.10526820438669089, -0.3423860463506388, -0.60158021270458772, -0.32020790415623457]),
        np.zeros((0,6)),
        np.array([]),
    )
    guess = np.array([0.4183733839886043, 1.685868545495764, 0.91029464030598595, 0.20609014773347709, 0.41772085105136075, 0.21627704946572152])
    cost, x = qps.find_local(guess, 100.0)
    qps.check_solution(x)


def test_local_case11():
    # The problem was that the last assert had a too low threhsold for the
    # current eps parameters.
    qps = QPSolver(
        np.array([
            [2.9180906444949501, -1.2533029767785395, 0.90478564555244811, 1.2787545043913675, 0.021113032381065766, -0.4116750982169069],
            [-1.2533029767785395, 5.6632299906311987, 1.4008760670203191, -2.0390934015965536, -0.85045793596990116, 0.73218455692376505],
            [0.90478564555244811, 1.4008760670203191, 5.2855381382892999, 2.677903087938867, 1.9188439544442164, -1.9709008178264937],
            [1.2787545043913675, -2.0390934015965536, 2.677903087938867, 4.993172076325008, 1.4299089422151234, -0.74354828220123381],
            [0.021113032381065766, -0.85045793596990116, 1.9188439544442164, 1.4299089422151234, 3.1450801244239486, -0.0036684167027846853],
            [-0.4116750982169069, 0.73218455692376505, -1.9709008178264937, -0.74354828220123381, -0.0036684167027846853, 2.5434283210738671],
        ]),
        np.array([-0.29192463967063931, -1.8024178730096532, 0.88758128602009223, 2.3598192316830855, -1.1550614364716303, 0.31271490996557316]),
        np.array([
            [-1.5989252225397905, 0.81244820675095797, -1.6307964553181504, -0.4456644455138748, 3.3167580671951278, -0.73081081535097092],
        ]),
        np.array([1.7184957352780219]),
        eps=1e-6,
    )
    cost0, x0 = qps.find_brute()
    qps.check_solution(x0)
    guess = np.array([0.5527936565693109, 1.2295737022405098, 0.68854432659855258, 0.11250382959013361, 0.84458517701828661, 0.0340219228211837])
    cost, x1 = qps.find_local(guess, 100.0)
    qps.check_solution(x1)
    assert abs(x0 - x1).max() < qps.eps
