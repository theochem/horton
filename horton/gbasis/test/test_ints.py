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

def test_fac2():
    assert fac2(-20) == 1
    assert fac2(0) == 1
    assert fac2(1) == 1
    assert fac2(2) == 2
    assert fac2(3) == 3
    assert fac2(4) == 8
    assert fac2(5) == 15


def test_binom():
    assert binom(1,1) == 1
    assert binom(5,3) == 10
    assert binom(3,2) == 3
    assert binom(10,4) == 210
    assert binom(18,14) == 3060
    assert binom(5, 1) == 5
    assert binom(5, 0) == 1
    assert binom(0, 0) == 1
    assert binom(5, 5) == 1


def test_gpt_coeff():
    def py_gpt_coeff(k, n0, n1, pa, pb):
        result = 0
        for q in xrange(max(-k, k-2*n0), min(k, 2*n1-k)+1, 2):
            i0 = (k+q)/2
            i1 = (k-q)/2
            assert (k+q)%2 == 0
            assert (k-q)%2 == 0
            result += binom(n0, i0)*binom(n1, i1)*pa**(n0-i0)*pb**(n1-i1)
        return result

    pa = 0.8769
    pb = 0.123
    for k in xrange(5):
        check = py_gpt_coeff(k, 2, 2, pa, pb)
        result = gpt_coeff(k, 2, 2, pa, pb)
        assert abs(check - result) < 1e-10

    for k in xrange(7):
        check = py_gpt_coeff(k, 3, 3, pa, pb)
        result = gpt_coeff(k, 3, 3, pa, pb)
        assert abs(check - result) < 1e-10


def test_gb_overlap_int1d():
    assert abs(gb_overlap_int1d(0, 0, 0.0, 0.0, 1/4.0) - 0.886227) < 1e-5
    print gb_overlap_int1d(2, 2, 0.0, 0.0, 1/5.0)
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/5.0) - 0.023780) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/5.0) - 0.023780) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.0, 0.0, 1/5.0) - 0.792665) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/1.0) - 1.329340) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/1.0) - 1.329340) < 1e-5
    assert abs(gb_overlap_int1d(1, 1, 0.0, 0.0, 1/1.0) - 0.886227) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/2.0) - 0.234996) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/2.0) - 0.234996) < 1e-5
    assert abs(gb_overlap_int1d(1, 1, 0.0, 0.0, 1/2.0) - 0.313329) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/3.0) - 0.085277) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/3.0) - 0.085277) < 1e-5
    assert abs(gb_overlap_int1d(1, 1, 0.0, 0.0, 1/3.0) - 0.170554) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/4.0) - 0.041542) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/4.0) - 0.041542) < 1e-5
    assert abs(gb_overlap_int1d(1, 1, 0.0, 0.0, 1/4.0) - 0.110778) < 1e-5
    assert abs(gb_overlap_int1d(2, 2, 0.0, 0.0, 1/5.0) - 0.023780) < 1e-5

    assert abs(gb_overlap_int1d(0, 0, 0.000000, -0.377945, 1/211400.020633) - 0.003855) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000000, -0.377945, 1/31660.020633) - 0.009961) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000001, -0.377944, 1/7202.020633) - 0.020886) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000004, -0.377941, 1/2040.020633) - 0.039243) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000012, -0.377934, 1/666.420633) - 0.068660) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000032, -0.377913, 1/242.020633) - 0.113933) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000082, -0.377864, 1/95.550633) - 0.181325) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000194, -0.377751, 1/40.250633) - 0.279376) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000440, -0.377506, 1/17.740633) - 0.420814) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, 0.000972, -0.376974, 1/8.025633) - 0.625656) < 1e-5

    assert abs(gb_overlap_int1d(0, 0, -0.014528, 0.363418, 1/8.325) - 0.614303) < 1e-5
    assert abs(gb_overlap_int1d(0, 3, 0.014528, -0.363418, 1/8.325) - -0.069710) < 1e-5
    assert abs(gb_overlap_int1d(0, 4, -0.014528, 0.363418, 1/8.325) - 0.046600) < 1e-5
    assert abs(gb_overlap_int1d(0, 3, -0.014528, 0.363418, 1/8.325) - 0.069710) < 1e-5
    assert abs(gb_overlap_int1d(0, 1, 0.014528, -0.363418, 1/8.325) - -0.223249) < 1e-5
    assert abs(gb_overlap_int1d(0, 0, -0.101693, 2.543923, 1/8.325) - 0.614303) < 1e-5
    assert abs(gb_overlap_int1d(0, 2, -0.014528, 0.363418, 1/8.325) - 0.118028) < 1e-5
    assert abs(gb_overlap_int1d(0, 1, -0.014528, 0.363418, 1/8.325) - 0.223249) < 1e-5
    assert abs(gb_overlap_int1d(0, 3, 0.014528, -0.363418, 1/8.325) - -0.069710) < 1e-5


def test_overlap_norm():
    max_nbasis = 15
    gb2oi = GB2OverlapIntegral(max_nbasis)
    for shell_type in 0, 1, 2, 3:
        for alpha in np.arange(0.5, 2.51, 0.5):
            scales = np.ones(get_shell_nbasis(shell_type), float)
            r = np.random.uniform(-1, 1, 3)
            gb2oi.reset(shell_type, shell_type, r, r)
            gb2oi.add(1.0, alpha, alpha, scales, scales)
            work = gb2oi.get_work()
            diag = np.diag(work)[:get_shell_nbasis(shell_type)]
            indexes = np.array([shell_type, 0, 0])
            counter = 0
            while True:
                check = diag[counter]*gob_normalization(alpha, indexes)**2
                assert abs(check - 1) < 1e-10
                if not iter_pow1_inc(indexes):
                    break
                counter += 1


def test_gb2_overlap_integral_class():
    max_nbasis = 15
    r0 = np.array([2.645617, 0.377945, -0.188973])
    r1 = np.array([0.000000, 0.000000, 0.188973])

    gb2oi = GB2OverlapIntegral(max_nbasis)
    assert gb2oi.max_nbasis == max_nbasis

    # case 1
    scales0 = np.ones(10, float)
    scales1 = np.ones(15, float)
    gb2oi.reset(3, 4, r0, r1)
    work = gb2oi.get_work()
    assert work.shape == (max_nbasis, max_nbasis)
    assert (work == 0.0).all()
    assert not (work is gb2oi.get_work())
    gb2oi.add(1.0, 5.398, 0.320, scales0, scales1)
    work = gb2oi.get_work()
    #assert abs(work[0, 3] - -0.0135976488493) < 1e-5
    assert abs(work[0, 3] - -0.001256) < 1e-5
    gb2oi.add(0.5, 5.398, 0.320, scales0, scales1)
    work = gb2oi.get_work()
    #assert abs(work[0, 3] - -1.5*0.0135976488493) < 1e-5
    assert abs(work[0, 3] - -1.5*0.001256) < 1e-5

    # case 2
    scales0 = np.ones(6, float)
    scales1 = np.ones(1, float)
    gb2oi.reset(2, 0, r0, r0)
    gb2oi.add(1.0, 0.463, 8.005, scales0, scales1)
    work = gb2oi.get_work()
    #assert abs(work[3,0] - 0.0193575315113) < 1e-5
    assert abs(work[3,0] - 0.013343) < 1e-5
    # 200 110 101 020 011 022
