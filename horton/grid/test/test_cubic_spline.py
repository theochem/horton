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


def test_tridiagsym_solve():
    N = 10
    A = np.zeros((N,N), float)
    # randomize the diagonal
    A.ravel()[::N+1] = np.random.uniform(1,2,N)
    # randomize the upper diagonal
    A.ravel()[1::N+1] = np.random.uniform(-1,0,N-1)
    # clone the lower diagonal
    A.ravel()[N::N+1] = A.ravel()[1::N+1]
    # test the inverse for all possible basis vectors
    Ainv = np.linalg.inv(A)
    for i in xrange(N):
        right = np.zeros(N, float)
        right[i] = 1.0
        solution = np.zeros(N, float)
        tridiagsym_solve(
            A.ravel()[::N+1].copy(),
            A.ravel()[N::N+1].copy(),
            right, solution
        )
        error = abs(solution - Ainv[:,i]).max()
        assert(error < 1e-12)


def test_tridiag_solve():
    N = 4
    A = np.zeros((N,N), float)
    # randomize the lower diagonal
    A.ravel()[N::N+1] = np.random.uniform(-1,1,N-1)
    # randomize the diagonal
    A.ravel()[::N+1] = np.random.uniform(1,2,N)
    # randomize the upper diagonal
    A.ravel()[1::N+1] = np.random.uniform(-1,1,N-1)
    # test the inverse for all possible basis vectors
    for i in xrange(N):
        right = np.zeros(N, float)
        right[i] = 1.0
        solution1 = np.zeros(N, float)
        tridiag_solve(
            A.ravel()[N::N+1].copy(),
            A.ravel()[::N+1].copy(),
            A.ravel()[1::N+1].copy(),
            right, solution1
        )
        right = np.zeros(N, float)
        right[i] = 1.0
        solution2 = np.linalg.solve(A, right)
        error = abs(solution1 - solution2).max()
        assert(error < 1e-12)


def test_basics():
    N = 10
    y = np.random.normal(0,1,N)
    d = np.random.normal(0,1,N)
    cs = CubicSpline(y,d)
    assert (cs.copy_y() == y).all()
    assert (cs.copy_d() == d).all()
    # test blindly function values at random grid points
    xnew = np.random.uniform(-5, 15, 50)
    ynew = np.zeros(len(xnew), float)
    cs(xnew, ynew)
    cs.integrate()


def test_continuity():
    N = 10
    y = np.random.normal(0,1,N)
    cs = CubicSpline(y)
    d = cs.copy_d()
    # test the function values at the grid points
    xnew = np.arange(N, dtype=float)
    ynew = np.zeros(len(xnew), float)
    cs(xnew, ynew)
    error = abs(y - ynew).max()
    assert(error < 1e-14)
    # test the second order derivative at the first point
    dd = 3*(y[1] - y[0]) - (d[1] + 2*d[0])
    error = abs(dd)
    assert(error < 1e-14)
    # test the second order derivative at the last point
    dd = 3*(y[N-2] - y[N-1]) + (d[N-2] + 2*d[N-1])
    error = abs(dd)
    assert(error < 1e-14)
    # test the continuity of the second order derivative
    for i in xrange(1,N-1):
        dd1 = 3*(y[i-1] - y[i]) + (d[i-1] + 2*d[i])
        dd2 = 3*(y[i+1] - y[i]) - (d[i+1] + 2*d[i])
        error = abs(dd1 - dd2)
        assert(error < 1e-14)


def test_accuracy():
    x_high = 2.0*np.pi
    size = 51
    delta = x_high/(size-1)
    t = np.arange(size, dtype=float)
    x = t*delta
    y = np.sin(x)
    cs = CubicSpline(y)
    newx = np.arange(0, 1.0005, 0.01)*x_high
    newt = newx/delta
    newy = np.zeros(len(newt))
    cs(newt, newy)
    error = abs(newy - np.sin(newx)).max()
    assert(error<1e-6)


def test_deriv_case1():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y, d)
    assert abs(cs.deriv(x) - d).max() < 1e-10
    x = np.arange(9, dtype=float)+0.5
    y = np.exp(-0.3*x)
    d = -0.3*y
    assert abs(cs.deriv(x) - d).max() < 1e-5


def test_deriv_case2():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y)
    assert abs(cs.deriv(x) - d).max() < 3e-2
    x = np.arange(9, dtype=float)+0.5
    y = np.exp(-0.3*x)
    d = -0.3*y
    assert abs(cs.deriv(x) - d).max() < 3e-2


def test_deriv_case3():
    y = np.random.normal(0, 1, 10)
    cs = CubicSpline(y)
    x = np.arange(9, dtype=float)+0.5
    eps = 1e-4
    d1 = cs.deriv(x)
    d2 = (cs(x+eps) - cs(x-eps))/(2*eps)
    assert abs(d1-d2).max() < 1e-6


def test_deriv_case4():
    y = np.random.normal(0, 1, 10)
    cs = CubicSpline(y)
    x = np.arange(9, dtype=float)+0.5
    eps = 1e-4
    d1 = cs.deriv(x)
    d2 = (cs(x+eps) - cs(x-eps))/(2*eps)
    assert abs(d1-d2).max() < 1e-6


def test_deriv2_case1():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y, d)
    d2 = 0.09*y
    assert abs(cs.deriv2(x) - d2).max() < 1e-3


def test_deriv2_case2():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y)
    assert abs(cs.deriv(x) - d).max() < 3e-2
    x = np.arange(9, dtype=float)+0.5
    y = np.exp(-0.3*x)
    d2 = 0.09*y
    assert abs(cs.deriv2(x) - d2)[2:-2].max() < 3e-3


def test_deriv2_case3():
    y = np.random.normal(0, 1, 10)
    cs = CubicSpline(y)
    x = np.arange(9, dtype=float)+0.5
    eps = 1e-4
    da = cs.deriv2(x)
    db = (cs.deriv(x+eps) - cs.deriv(x-eps))/(2*eps)
    assert abs(da-db).max() < 1e-6


def test_deriv2_case4():
    y = np.random.normal(0, 1, 10)
    cs = CubicSpline(y)
    x = np.arange(9, dtype=float)+0.5
    eps = 1e-4
    da = cs.deriv2(x)
    db = (cs.deriv(x+eps) - cs.deriv(x-eps))/(2*eps)
    assert abs(da-db).max() < 1e-6


def test_extrapolation():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y, d)
    newx = np.array([-2.5, -1.1, 10.5, 11.5])
    assert abs(cs(newx) - np.exp(-0.3*newx)).max() < 1e-10
