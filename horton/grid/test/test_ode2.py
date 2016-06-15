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


from scipy.special import erf
import numpy as np, random

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


def _merge(y, d):
    '''Put y and d in one vector'''
    return np.array([y, d]).T.ravel()


def test_hermite_overlap2():
    # things that should be zero (no overlap)
    assert hermite_overlap2(10, 2*5+0, False, 2*2+0, False) == 0
    assert hermite_overlap2(10, 2*2+0, False, 2*5+0, False) == 0
    # things that should be zero (uneven integrand)
    assert hermite_overlap2(10, 2*5+0, False, 2*5+1, False) == 0
    assert hermite_overlap2(10, 2*5+1, False, 2*5+0, False) == 0
    assert hermite_overlap2(10, 2*5+0, False, 2*5+0, True ) == 0
    assert hermite_overlap2(10, 2*5+0, True , 2*5+0, False) == 0
    # things that should be non-zero
    assert hermite_overlap2(10, 2*5+0, False, 2*5+0, False) != 0
    assert hermite_overlap2(10, 2*5+0, True , 2*5+0, True ) != 0
    assert hermite_overlap2(10, 2*5+1, False, 2*5+1, False) != 0
    assert hermite_overlap2(10, 2*5+1, True , 2*5+1, True ) != 0
    # symmetry
    assert hermite_overlap2(10, 2*5+0, False, 2*6+0, False) == hermite_overlap2(10, 2*6+0, False, 2*5+0, False)
    assert hermite_overlap2(10, 2*5+0, True , 2*6+0, True ) == hermite_overlap2(10, 2*6+0, True , 2*5+0, True )
    # anti-symmetry
    assert hermite_overlap2(10, 2*5+0, False, 2*6+1, False) == -hermite_overlap2(10, 2*6+0, False, 2*5+1, False)
    assert hermite_overlap2(10, 2*5+0, True , 2*6+1, True ) == -hermite_overlap2(10, 2*6+0, True , 2*5+1, True )
    # xmax
    assert 2*hermite_overlap2(8, 2*0+0, False, 2*0+0, False) == hermite_overlap2(8, 2*4+0, False, 2*4+0, False)
    assert 2*hermite_overlap2(8, 2*8+0, False, 2*8+0, False) == hermite_overlap2(8, 2*4+0, False, 2*4+0, False)


def test_hermite_overlap3():
    # things that should be zero
    assert hermite_overlap3(10, 2*5+0, False, 2*5+0, False, 2*7+0, False) == 0
    assert hermite_overlap3(10, 2*5+0, False, 2*6+0, False, 2*7+0, False) == 0
    # things that should be zero (uneven integrand)
    assert hermite_overlap3(10, 2*5+0, False, 2*5+0, False, 2*5+0, True ) == 0
    assert hermite_overlap3(10, 2*5+0, False, 2*5+1, False, 2*5+0, False) == 0
    # things that should be non-zero (even integrand)
    assert hermite_overlap3(10, 2*5+0, False, 2*5+0, False, 2*5+0, False) != 0
    assert hermite_overlap3(10, 2*5+0, True , 2*5+0, False, 2*5+0, True ) != 0
    assert hermite_overlap3(10, 2*5+1, False, 2*5+0, False, 2*5+1, False) != 0


def test_hermite_node():
    # key function values and derivatives for kind=0
    assert hermite_node(0, 0, 0, False) == 1
    assert hermite_node(1, 0, 0, False) == 0
    assert hermite_node(-1, 0, 0, False) == 0
    assert hermite_node(2, 0, 0, False) == 0
    assert hermite_node(-2, 0, 0, False) == 0
    assert hermite_node(0, 0, 0, True) == 0
    assert hermite_node(1, 0, 0, True) == 0
    assert hermite_node(-1, 0, 0, True) == 0
    assert hermite_node(2, 0, 0, True) == 0
    assert hermite_node(-2, 0, 0, True) == 0
    # key function values and derivatives for kind=1
    assert hermite_node(0, 0, 1, False) == 0
    assert hermite_node(1, 0, 1, False) == 0
    assert hermite_node(-1, 0, 1, False) == 0
    assert hermite_node(2, 0, 1, False) == 0
    assert hermite_node(-2, 0, 1, False) == 0
    assert hermite_node(0, 0, 1, True) == 1
    assert hermite_node(1, 0, 1, True) == 0
    assert hermite_node(-1, 0, 1, True) == 0
    assert hermite_node(2, 0, 1, True) == 0
    assert hermite_node(-2, 0, 1, True) == 0
    # outside the interval
    for kind in 0, 1:
        for deriv in False, True:
            assert hermite_node(2, 0, kind, deriv) == 0
            assert hermite_node(-2, 0, kind, deriv) == 0
            assert hermite_node(8, 10, kind, deriv) == 0
            assert hermite_node(12, 10, kind, deriv) == 0


def test_hermite_product2():
    assert hermite_product2(0, 2*0+0, False, 2*0+0, False) == 1
    assert hermite_product2(0, 2*1+0, False, 2*0+0, False) == 0
    assert hermite_product2(0, 2*0+0, False, 2*0+1, False) == 0


def test_build_ode2_basics():
    npoint = 4
    by = np.random.normal(0, 1, npoint)
    bd = np.random.normal(0, 1, npoint)
    ay = np.random.normal(0, 1, npoint)
    ad = np.random.normal(0, 1, npoint)
    fy = np.random.normal(0, 1, npoint)
    fd = np.random.normal(0, 1, npoint)
    bcs = [None, None, np.random.normal(0, 1), np.random.normal(0, 1)]
    random.shuffle(bcs)
    coeffs, rhs = build_ode2(by, bd, ay, ad, fy, fd, bcs)
    assert coeffs.shape == (2*npoint, 2*npoint)
    for irow, icol in (0, 4), (0, 5), (0, 6), (0, 7), (2, 6), (2, 7):
        assert coeffs[irow, icol] == 0
        assert coeffs[irow+1, icol] == 0
        assert coeffs[icol, irow] == 0
        assert coeffs[icol, irow+1] == 0
    assert rhs.shape == (2*npoint,)
    assert abs(rhs).min() != 0


def test_simple_quadratic():
    npoint = 2
    by = np.zeros(npoint)
    bd = np.zeros(npoint)
    ay = np.zeros(npoint)
    ad = np.zeros(npoint)
    fy = np.zeros(npoint) + 1
    fd = np.zeros(npoint)
    coeffs, rhs = build_ode2(by, bd, ay, ad, fy, fd, (0.0, None, 0.0, None))
    # case 1
    solution = np.linalg.solve(coeffs, rhs)
    assert abs(solution - [0.0, -0.5, 0.0, 0.5]).max() < 1e-12
    # case 2
    rhs[0] = 1.0
    rhs[-2] = 1.0
    solution = np.linalg.solve(coeffs, rhs)
    assert abs(solution - [1.0, -0.5, 1.0, 0.5]).max() < 1e-12
    # case 3
    rhs[0] = 0.0
    rhs[-2] = 0.5
    solution = np.linalg.solve(coeffs, rhs)
    assert abs(solution - [0.0, 0.0, 0.5, 1.0]).max() < 1e-12
    # case 4
    rhs[0] = 2.0
    rhs[-2] = 0.5
    solution = np.linalg.solve(coeffs, rhs)
    assert abs(solution - [2.0, -2.0, 0.5, -1.0]).max() < 1e-12


def test_long_quadratic():
    npoint = 4
    by = np.zeros(npoint)
    bd = np.zeros(npoint)
    ay = np.zeros(npoint)
    ad = np.zeros(npoint)
    fy = np.zeros(npoint) + 1
    fd = np.zeros(npoint)
    x = np.arange(npoint, dtype=float)
    coeffs, rhs = build_ode2(by, bd, ay, ad, fy, fd, (0.0, None, 0.0, None))
    # case 1
    solution = np.linalg.solve(coeffs, rhs)
    y = 0.5*(x**2 - 3.0*x)
    d = x - 1.5
    assert abs(solution - _merge(y, d)).max() < 1e-12
    # case 2
    rhs[0] = 1.0
    rhs[-2] = 1.0
    solution = np.linalg.solve(coeffs, rhs)
    y = 0.5*(x**2 - 3.0*x) + 1.0
    d = x - 1.5
    assert abs(solution - _merge(y, d)).max() < 1e-12
    # case 3
    rhs[0] = 0.0
    rhs[-2] = 4.5
    solution = np.linalg.solve(coeffs, rhs)
    y = 0.5*x**2
    d = x
    assert abs(solution - _merge(y, d)).max() < 1e-12
    # case 4
    rhs[0] = 2.0
    rhs[-2] = 0.5
    solution = np.linalg.solve(coeffs, rhs)
    y = 0.5*x**2 - 2*x + 2
    d = x - 2
    assert abs(solution - _merge(y, d)).max() < 1e-12


def test_simple_cubic():
    npoint = 2
    by = np.zeros(npoint)
    bd = np.zeros(npoint)
    ay = np.zeros(npoint)
    ad = np.zeros(npoint)
    fy = np.array([0.0, 1.0])
    fd = np.array([1.0, 1.0])
    coeffs, rhs = build_ode2(by, bd, ay, ad, fy, fd, (0.0, None, 0.0, None))
    # case 1
    solution = np.linalg.solve(coeffs, rhs)
    assert abs(solution - [0, -1.0/6.0, 0, 1.0/3.0]).max() < 1e-12
    # case 2
    rhs[0] = 1.0
    rhs[-2] = 1.0
    solution = np.linalg.solve(coeffs, rhs)
    assert abs(solution - [1, -1.0/6.0, 1, 1.0/3.0]).max() < 1e-12
    # case 3
    rhs[0] = 0.0
    rhs[-2] = 1.0/6.0
    solution = np.linalg.solve(coeffs, rhs)
    assert abs(solution - [0.0, 0.0, 1.0/6.0, 0.5]).max() < 1e-12


def test_long_cubic():
    npoint = 4
    x = np.arange(npoint, dtype=float)
    by = np.zeros(npoint)
    bd = np.zeros(npoint)
    ay = np.zeros(npoint)
    ad = np.zeros(npoint)
    fy = x
    fd = np.ones(npoint)
    coeffs, rhs = build_ode2(by, bd, ay, ad, fy, fd, (0.0, None, 0.0, None))
    # case 1
    solution = np.linalg.solve(coeffs, rhs)
    y = x**3/6.0 - 1.5*x
    d = x**2/2.0 - 1.5
    assert abs(solution - _merge(y, d)).max() < 1e-12
    # case 2
    rhs[0] = 1.0
    rhs[-2] = 1.0
    solution = np.linalg.solve(coeffs, rhs)
    y = x**3/6.0 - 1.5*x + 1.0
    d = x**2/2.0 - 1.5
    assert abs(solution - _merge(y, d)).max() < 1e-12
    # case 3
    rhs[0] = 0.0
    rhs[-2] = 4.5
    solution = np.linalg.solve(coeffs, rhs)
    y = x**3/6.0
    d = x**2/2.0
    assert abs(solution - _merge(y, d)).max() < 1e-12
    # case 4
    rhs[0] = 2.0
    rhs[-2] = 0.5
    solution = np.linalg.solve(coeffs, rhs)
    y = x**3/6.0 - 2.0*x + 2.0
    d = x**2/2.0 - 2.0
    assert abs(solution - _merge(y, d)).max() < 1e-12


def test_poisson_s_identity():
    npoint = 50
    x = np.arange(0, npoint, dtype=float)
    sigma = 15.0
    rhoy = np.exp(-0.5*(x/sigma)**2)/sigma**3/(2*np.pi)**1.5
    rhod = np.exp(-0.5*(x/sigma)**2)/sigma**3/(2*np.pi)**1.5*(-x/sigma)/sigma

    # input for the solver
    fy = -4*np.pi*x*rhoy
    fd = -4*np.pi*(rhoy + x*rhod)
    by = np.zeros(npoint)
    bd = np.zeros(npoint)
    ay = np.zeros(npoint)
    ad = np.zeros(npoint)
    #ay = -1.0/x**2
    #ad = 2.0/x**3

    # get matrices
    coeffs, rhs = build_ode2(by, bd, ay, ad, fy, fd, (0.0, None, 1.0, None))
    solution = np.linalg.solve(coeffs, rhs)

    soly = erf(x/(np.sqrt(2)*sigma))
    sold = np.exp(-(x/(np.sqrt(2)*sigma))**2)*2/np.sqrt(np.pi)/(np.sqrt(2)*sigma)

    assert abs(solution[::2] - soly).max() < 2e-3
    assert abs(solution[1::2] - sold).max() < 2e-3


def test_poisson_s_exp():
    npoint = 50
    rtf = ExpRTransform(1e-3, 1e1, npoint)
    x = rtf.get_radii()
    sigma = 1.0
    rhoy = np.exp(-0.5*(x/sigma)**2)/sigma**3/(2*np.pi)**1.5
    rhod = np.exp(-0.5*(x/sigma)**2)/sigma**3/(2*np.pi)**1.5*(-x/sigma)/sigma

    # input for the solver
    fy = -4*np.pi*x*rhoy
    fd = -4*np.pi*(rhoy + x*rhod)

    # transform to linear coordinate
    j1 = rtf.get_deriv()
    j2 = rtf.get_deriv2()
    j3 = rtf.get_deriv3()
    by1 = -j2/j1
    bd1 = (j2*j2 - j1*j3)/(j1*j1)
    ay1 = np.zeros(npoint)
    ad1 = np.zeros(npoint)
    fy1 = fy*j1*j1
    fd1 = (fd*j1*j1 + 2*j2*fy)*j1

    # get matrices
    coeffs, rhs = build_ode2(by1, bd1, ay1, ad1, fy1, fd1, (0.0, None, 1.0, None))
    solution = np.linalg.solve(coeffs, rhs)
    solution[1::2] /= j1

    soly = erf(x/(np.sqrt(2)*sigma))
    sold = np.exp(-(x/(np.sqrt(2)*sigma))**2)*2/np.sqrt(np.pi)/(np.sqrt(2)*sigma)

    assert abs(solution[::2] - soly).max() < 2e-3
    assert abs(solution[1::2] - sold).max() < 2e-3


def check_solve_s(rtf, sigma=1.0, epsy=1e-6, epsd=1e-3):
    x = rtf.get_radii()
    z = np.zeros(len(x))

    b = CubicSpline(z, z, rtf)
    a = CubicSpline(z, z, rtf)

    rhoy = np.exp(-0.5*(x/sigma)**2)/sigma**3/(2*np.pi)**1.5
    rhod = np.exp(-0.5*(x/sigma)**2)/sigma**3/(2*np.pi)**1.5*(-x/sigma)/sigma
    fy = -4*np.pi*x*rhoy
    fd = -4*np.pi*(rhoy + x*rhod)
    f = CubicSpline(fy, fd, rtf)

    bcs = (0.0, None, 1.0, None)

    u = solve_ode2(b, a, f, bcs)

    soly = erf(x/(np.sqrt(2)*sigma))
    sold = np.exp(-(x/(np.sqrt(2)*sigma))**2)*2/np.sqrt(np.pi)/(np.sqrt(2)*sigma)
    assert abs(u.y - soly).max()/abs(soly).max() < epsy
    assert abs(u.dx - sold).max()/abs(sold).max() < epsd


def test_solve_ode2_s_identity():
    check_solve_s(IdentityRTransform(50), 10.0)


def test_solve_ode2_s_linear():
    check_solve_s(LinearRTransform(0.0, 4e1, 50), 8.0)


def test_solve_ode2_s_exp():
    check_solve_s(ExpRTransform(1e-6, 4e1, 100), 8.0, 1e-5, 1e-4)


def test_solve_ode2_s_power():
    check_solve_s(PowerRTransform(4e-4, 4e1, 100), 8.0, 1e-4, 1e-2)


def check_solve_s2(rtf, sigma=1.0, epsy=1e-10, epsd=1e-10):
    r = rtf.get_radii()
    z = np.zeros(len(r))

    b = CubicSpline(2/r, -2/r**2, rtf)
    a = CubicSpline(z, z, rtf)

    rhoy = np.exp(-0.5*(r/sigma)**2)/sigma**3/(2*np.pi)**1.5
    rhod = np.exp(-0.5*(r/sigma)**2)/sigma**3/(2*np.pi)**1.5*(-r/sigma)/sigma
    fy = -4*np.pi*rhoy
    fd = -4*np.pi*rhod
    f = CubicSpline(fy, fd, rtf)

    bcs = (None, 0.0, 1.0/r[-1], None)

    u = solve_ode2(b, a, f, bcs)


    s2s = np.sqrt(2)*sigma
    soly = erf(r/s2s)/r
    sold = np.exp(-(r/s2s)**2)*2/np.sqrt(np.pi)/s2s/r - erf(r/s2s)/r**2
    if False:
        print abs(u.y - soly).max()/abs(soly).max()
        print abs(u.dx - sold).max()/abs(sold).max()
        import matplotlib.pyplot as pt
        pt.clf()
        pt.plot(r, u.y, label='numerical')
        pt.plot(r, soly, label='exact')
        pt.legend(loc=0)
        pt.savefig('s2.png')
    assert abs(u.y - soly).max()/abs(soly).max() < epsy
    assert abs(u.dx - sold).max()/abs(sold).max() < epsd


def test_solve_ode2_s2_linear():
    check_solve_s2(LinearRTransform(1e-5, 4e1, 200), 8.0, 1e-3, 4e-2)


def test_solve_ode2_s2_exp():
    check_solve_s2(ExpRTransform(1e-6, 4e1, 100), 8.0, 1e-5, 1e-3)


def test_solve_ode2_s2_power():
    check_solve_s2(PowerRTransform(4e-4, 4e1, 100), 8.0, 1e-6, 1e-4)


def check_solve_xexp(rtf, s=1.0):
    x = rtf.get_radii()
    o = np.ones(len(x))

    b = CubicSpline(s*(2 + s*x), o*s**2, rtf)
    a = CubicSpline(-s**3*x, -o*s**3, rtf)
    f = CubicSpline(4*s*(1+s*x)*np.exp(s*x), 4*s**2*(2+s*x)*np.exp(s*x), rtf)

    bcs_base = [
        x[0]*np.exp(s*x[0]),
        (1+s*x[0])*np.exp(s*x[0]),
        x[-1]*np.exp(s*x[-1]),
        (1+s*x[-1])*np.exp(s*x[-1]),
    ]
    for i0 in xrange(2):
        for i1 in xrange(2, 4):
            bcs = list(bcs_base)
            bcs[i0] = None
            bcs[i1] = None

            u = solve_ode2(b, a, f, bcs)

            soly = x*np.exp(s*x)
            sold = (1+s*x)*np.exp(s*x)
            if False:
                print i0, i1,
                print abs(u.y - soly).max()/abs(soly).max(),
                print abs(u.dx - sold).max()/abs(sold).max()
                print '%+10.5f  %+10.5f  %s' % (u.y[0], soly[0], bcs[0])
                print '%+10.5f  %+10.5f  %s' % (u.dx[0], sold[0], bcs[1])
                print '%+10.5f  %+10.5f  %s' % (u.y[-1], soly[-1], bcs[2])
                print '%+10.5f  %+10.5f  %s' % (u.dx[-1], sold[-1], bcs[3])
                import matplotlib.pyplot as pt
                pt.clf()
                pt.plot(x, u.y, label='numerical')
                pt.plot(x, soly, label='exact')
                pt.legend(loc=0)
                pt.savefig('xexp_%i_%i.png' % (i0, i1))
            assert abs(u.y - soly).max()/abs(soly).max() < 1e-3
            assert abs(u.dx - sold).max()/abs(sold).max() < 1e-3


def test_solve_ode2_xexp_identity():
    check_solve_xexp(IdentityRTransform(50), 0.01)


def test_solve_ode2_xexp_linear():
    check_solve_xexp(LinearRTransform(0.2, 1.2, 50))


def test_solve_ode2_xexp_exp():
    check_solve_xexp(ExpRTransform(0.2, 1.2, 50))


def test_solve_ode2_xexp_power():
    check_solve_xexp(PowerRTransform(0.0002, 2.0, 50))
