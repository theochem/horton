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
import h5py as h5

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import

from horton.test.common import numpy_seed


def test_tridiagsym_solve():
    N = 10
    A = np.zeros((N,N), float)
    for irep in xrange(100):
        with numpy_seed(irep):  # Make random numbers reproducible
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
            assert(error < 1e-9)


def test_basics_identity():
    N = 10
    y = np.random.normal(0,1,N)
    d = np.random.normal(0,1,N)
    cs = CubicSpline(y,d)
    assert (cs.y == y).all()
    assert (cs.dx == d).all()
    assert (cs.dt == d).all()


def test_basics_linear():
    N = 10
    y = np.random.normal(0,1,N)
    d = np.random.normal(0,1,N)
    rtf = LinearRTransform(-0.3, 0.6, N)
    cs = CubicSpline(y, d, rtf)
    assert (cs.y == y).all()
    assert (cs.dx == d).all()
    dt = d*rtf.get_deriv()
    assert abs(rtf.get_deriv() - 0.1).max() < 1e-10
    assert abs(cs.dt - dt).max() < 1e-15


def test_basics_exp():
    N = 10
    y = np.random.normal(0,1,N)
    d = np.random.normal(0,1,N)
    rtf = ExpRTransform(0.1, 1.0, N)
    cs = CubicSpline(y, d, rtf)
    assert (cs.y == y).all()
    assert (cs.dx == d).all()
    dt = d*rtf.get_deriv()
    assert abs(cs.dt - dt).max() < 1e-15


def check_continuity(ynew, y, d, N):
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


def test_continuity_identity():
    N = 10
    y = np.random.normal(0,1,N)
    cs = CubicSpline(y)
    # test the function values at the grid points
    xnew = np.arange(N, dtype=float)
    ynew = np.zeros(len(xnew), float)
    cs(xnew, ynew)
    check_continuity(ynew, y, cs.dt, N)


def test_continuity_exp():
    N = 10
    rtf = ExpRTransform(0.1, 1.0, N)
    y = np.random.normal(0,1,N)
    cs = CubicSpline(y, rtransform=rtf)
    # test the function values at the grid points
    tnew = np.arange(N, dtype=float)
    xnew = rtf.radius(tnew)
    ynew = np.zeros(len(xnew), float)
    cs(xnew, ynew)
    check_continuity(ynew, y, cs.dt, N)


def test_accuracy_identity():
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


def test_accuracy_exp():
    size = 51
    rtf = ExpRTransform(0.1, 1.0, size)
    x = rtf.get_radii()
    y = np.sin(x)
    cs = CubicSpline(y, rtransform=rtf)
    newx = np.arange(0.1, 1.0, 0.01)
    newy = np.zeros(len(newx))
    cs(newx, newy)
    error = abs(newy - np.sin(newx)).max()
    assert(error<4e-5)


def test_deriv_identity1():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y, d)
    assert abs(cs.deriv(x) - d).max() < 1e-10
    x = np.arange(9, dtype=float)+0.5
    y = np.exp(-0.3*x)
    d = -0.3*y
    assert abs(cs.deriv(x) - d).max() < 1e-5


def test_deriv_exp1():
    rtf = ExpRTransform(0.1, 1.0, 10)
    x = rtf.get_radii()
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y, d, rtf)
    assert abs(cs.deriv(x) - d).max() < 1e-15
    t = np.arange(9, dtype=float)+0.5
    x = rtf.radius(t)
    y = np.exp(-0.3*x)
    d = -0.3*y
    assert abs(cs.deriv(x) - d).max() < 1e-6


def test_deriv_identity2():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y)
    assert abs(cs.deriv(x) - d).max() < 3e-2
    x = np.arange(9, dtype=float)+0.5
    y = np.exp(-0.3*x)
    d = -0.3*y
    assert abs(cs.deriv(x) - d).max() < 3e-2


def test_deriv_exp2():
    rtf = ExpRTransform(0.1, 1.0, 10)
    x = rtf.get_radii()
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y, rtransform=rtf)
    assert abs(cs.deriv(x) - d).max() < 3e-2
    t = np.arange(9, dtype=float)+0.5
    x = rtf.radius(t)
    y = np.exp(-0.3*x)
    d = -0.3*y
    assert abs(cs.deriv(x) - d).max() < 4e-3


def test_deriv_identity3():
    y = np.random.normal(0, 1, 10)
    cs = CubicSpline(y)
    x = np.arange(9, dtype=float)+0.5
    eps = 1e-4
    d1 = cs.deriv(x)
    d2 = (cs(x+eps) - cs(x-eps))/(2*eps)
    assert abs(d1-d2).max() < 1e-6


def test_deriv_exp3():
    rtf = ExpRTransform(0.1, 1.0, 10)
    y = np.random.normal(0, 1, 10)
    cs = CubicSpline(y, rtransform=rtf)
    t = np.arange(9, dtype=float)+0.5
    x = rtf.radius(t)
    eps = 1e-6
    d1 = cs.deriv(x)
    d2 = (cs(x+eps) - cs(x-eps))/(2*eps)
    assert abs(d1-d2).max() < 1e-6


def test_deriv_identity4():
    y = np.random.normal(0, 1, 10)
    d = np.random.normal(0, 1, 10)
    cs = CubicSpline(y, d)
    x = np.arange(9, dtype=float)+0.5
    eps = 1e-4
    d1 = cs.deriv(x)
    d2 = (cs(x+eps) - cs(x-eps))/(2*eps)
    assert abs(d1-d2).max() < 1e-6


def test_deriv_exp4():
    rtf = ExpRTransform(0.1, 1.0, 10)
    y = np.random.normal(0, 1, 10)
    d = np.random.normal(0, 1, 10)
    cs = CubicSpline(y, d, rtf)
    t = np.arange(9, dtype=float)+0.5
    x = rtf.radius(t)
    eps = 1e-6
    d1 = cs.deriv(x)
    d2 = (cs(x+eps) - cs(x-eps))/(2*eps)
    assert abs(d1-d2).max() < 1e-6


def test_extrapolation1_identity():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y, d)
    newx = np.array([-2.5, -1.1])
    assert abs(cs(newx) - np.exp(-0.3*newx)).max() < 1e-10
    assert abs(cs.deriv(newx) - -0.3*np.exp(-0.3*newx)).max() < 1e-10
    newx = np.array([10.5, 11.5])
    assert abs(cs(newx)).max() < 1e-10
    assert abs(cs.deriv(newx)).max() < 1e-10


def test_extrapolation2_identity():
    x = np.arange(10, dtype=float)
    y = x**2 + 1
    d = 2*x
    cs = CubicSpline(y, d)
    newx = np.array([-2.5, -1.1])
    assert abs(cs(newx) - 1.0).max() < 1e-10
    assert abs(cs.deriv(newx)).max() < 1e-10


def test_extrapolation1_exp():
    rtf = ExpRTransform(0.1, 1.0, 10)
    x = rtf.get_radii()
    y = np.exp(-0.3*x)
    d = -0.3*y
    cs = CubicSpline(y, d, rtf)
    newx = np.array([0.001, 0.01])
    assert abs(cs(newx) - np.exp(-0.3*newx)).max() < 1e-10
    assert abs(cs.extrapolation.eval_left(newx[0]) - np.exp(-0.3*newx[0])).max() < 1e-10
    assert abs(cs.deriv(newx) - -0.3*np.exp(-0.3*newx)).max() < 1e-10
    assert abs(cs.extrapolation.deriv_left(newx[0]) - -0.3*np.exp(-0.3*newx[0])).max() < 1e-10
    newx = np.array([1.1, 10.1])
    assert abs(cs(newx)).max() < 1e-10
    assert abs(cs.deriv(newx)).max() < 1e-10


def test_bug_steven():
    y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.95307674e-005, 2.64520014e-004, 4.86536423e-004, 7.72847240e-004, 1.15087885e-003, 1.66508373e-003, 2.38709713e-003, 3.42860950e-003, 4.95314131e-003, 7.17684800e-003, 1.03407874e-002, 1.46341247e-002, 2.00627694e-002, 2.63026157e-002, 3.26349709e-002, 3.80678816e-002, 4.16324963e-002, 4.26735205e-002, 4.09493929e-002, 3.65937608e-002, 3.01384965e-002, 2.25802551e-002, 1.52220562e-002, 9.20415761e-003, 5.04355509e-003, 2.56723646e-003, 1.24504774e-003, 5.72945090e-004, 2.36808736e-004, 8.00848384e-005, 1.99672917e-005, 3.30308018e-006, 3.22840493e-007, 1.62065171e-008, 3.50644056e-010, 2.62139812e-012, 5.12409163e-015, 1.84238101e-018, 7.81842741e-023, 2.23721221e-028, 2.12999970e-035, 2.76742060e-044, 1.59368870e-055, 9.84245631e-070, 1.08780738e-087, 2.24648735e-110])
    rtf = ExpRTransform(0.0003779452267842504, 37.79452267842504, 100)
    cs = CubicSpline(y, rtransform=rtf)
    d = np.array([0.000141, 0.00141, 0.0141, 0.141, 1.41])
    s = cs(d)
    assert not np.isnan(s).any()


def test_extrapolation_identity_power():
    x = np.arange(10, dtype=float)
    y = 1/(np.exp(-2*x)+x**2)
    d = -y**2*(np.exp(-2*x)*(-2) + 2*x)
    cs = CubicSpline(y, d, extrapolation=PowerExtrapolation(-2.0))
    assert cs.extrapolation.power == -2.0
    newx = np.array([-2.5, -1.1])
    assert abs(cs(newx)).max() == 0.0
    assert abs(cs.deriv(newx)).max() == 0.0
    newx = np.array([10.5, 11.5])
    assert abs(cs(newx) - 1/newx**2).max() < 1e-10
    assert abs(cs.extrapolation.eval_right(newx[0]) - 1/newx[0]**2).max() < 1e-10
    assert abs(cs.deriv(newx) - -2/newx**3).max() < 1e-10
    assert abs(cs.extrapolation.deriv_right(newx[0]) - -2/newx[0]**3) < 1e-10


def test_extrapolation_exp_power():
    rtf = ExpRTransform(0.1, 1.0, 10)
    x = rtf.get_radii()
    y = x**2
    d = 2*x
    cs = CubicSpline(y, d, rtf, PowerExtrapolation(2))
    assert cs.extrapolation.power == 2.0
    newx = np.array([0.001, 0.01])
    assert abs(cs(newx)).max() == 0.0
    assert abs(cs.deriv(newx)).max() == 0.0
    newx = np.array([1.1, 10.1])
    assert abs(cs(newx) - newx**2).max() < 1e-10
    assert abs(cs.extrapolation.eval_right(newx[0]) - newx[0]**2).max() < 1e-10
    assert abs(cs.deriv(newx) - 2*newx).max() < 1e-10
    assert abs(cs.extrapolation.deriv_right(newx[0]) - 2*newx[0]) < 1e-10


def test_extrapolation_exp_potential():
    rtf = ExpRTransform(0.1, 1.0, 10)
    x = rtf.get_radii()
    y = 1/(np.exp(-2*x)+x**2)
    d = -y**2*(np.exp(-2*x)*(-2) + 2*x)
    for l in 0, 1, 2, 3:
        cs = CubicSpline(y, d, rtf, extrapolation=PotentialExtrapolation(l))
        assert cs.extrapolation.l == l
        np.testing.assert_allclose(cs.extrapolation.amp_left, y[0]/x[0]**l)
        np.testing.assert_allclose(cs.extrapolation.amp_right, y[-1]*x[-1]**(l+1))
        newx = np.array([-2.5, -1.1])
        amp_left = cs.extrapolation.amp_left
        np.testing.assert_allclose(cs(newx), amp_left*newx**l)
        np.testing.assert_allclose(cs.deriv(newx), amp_left*l*newx**(l-1))
        newx = np.array([10.5, 11.5])
        amp_right = cs.extrapolation.amp_right
        np.testing.assert_allclose(cs(newx), amp_right/newx**(l+1))
        np.testing.assert_allclose(cs.deriv(newx), -(l+1)*amp_right/newx**(l+2))


def test_consistency_h5():
    with h5.File('horton.grid.test.test_cubic_spline.test_consistency_h5', driver='core', backing_store=False) as chk:
        rtf = ExpRTransform(0.1, 1.0, 10)
        y = np.random.normal(0, 1, 10)
        d = np.random.normal(0, 1, 10)
        cs1 = CubicSpline(y, d, rtf, PowerExtrapolation(2))
        cs1.to_hdf5(chk)
        cs2 = CubicSpline.from_hdf5(chk)
        assert (cs1.y == cs2.y).all()
        assert (cs1.dx == cs2.dx).all()
        assert (cs1.dt == cs2.dt).all()
        assert cs1.rtransform.to_string() == cs2.rtransform.to_string()
        assert cs1.extrapolation.to_string() == cs2.extrapolation.to_string()


def test_zero_extrapolation_string():
    assert ZeroExtrapolation().to_string() == 'ZeroExtrapolation'
    assert isinstance(Extrapolation.from_string('ZeroExtrapolation'), ZeroExtrapolation)


def test_cusp_extrapolation_string():
    assert CuspExtrapolation().to_string() == 'CuspExtrapolation'
    assert isinstance(Extrapolation.from_string('CuspExtrapolation'), CuspExtrapolation)


def test_power_extrapolation_string():
    assert PowerExtrapolation(2.0).to_string() == 'PowerExtrapolation 2.0'
    ep = Extrapolation.from_string(PowerExtrapolation(5.1247953315476).to_string())
    assert isinstance(ep, PowerExtrapolation)
    assert ep.power == 5.1247953315476
