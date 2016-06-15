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


try:
    from sympy import S, sqrt
except ImportError:
    from nose.plugins.skip import SkipTest
    raise SkipTest

from harmonics import get_poly_conversion


# Comparison with some results from:
# Schlegel, H. B.;  Frisch, M. J. Int. J. Quantum Chem. 1995, 54, 83-87.

def test_lcs_d():
    lcs2 = get_poly_conversion(2)
    assert lcs2[0,0] == -1/S(2) # xx
    assert lcs2[0,1] == 0 # xy
    assert lcs2[0,2] == 0 # xz
    assert lcs2[0,3] == -1/S(2) # yy
    assert lcs2[0,4] == 0 # yz
    assert lcs2[0,5] == 1 # zz

def test_lcs_f():
    lcs3 = get_poly_conversion(3)
    assert lcs3[0,0] == 0 # xxx
    assert lcs3[0,1] == 0 # xxy
    assert lcs3[0,2] == -3/S(2)/sqrt(5) # xxz
    assert lcs3[0,3] == 0 # xyy
    assert lcs3[0,4] == 0 # xyz
    assert lcs3[0,5] == 0 # xzz
    assert lcs3[0,6] == 0 # yyy
    assert lcs3[0,7] == -3/S(2)/sqrt(5) # yyz
    assert lcs3[0,8] == 0 # yzz
    assert lcs3[0,9] == 1 # zzz

def test_lcs_g():
    lcs4 = get_poly_conversion(4)
    assert lcs4[0,0] == 3/S(8) # xxxx
    assert lcs4[0,1] == 0 # xxxy
    assert lcs4[0,2] == 0 # xxxz
    assert lcs4[0,3] == 3*sqrt(3)/sqrt(35)/S(4) # xxyy
    assert lcs4[0,4] == 0 # xxyz
    assert lcs4[0,5] == -3*sqrt(3)/sqrt(35) # xxzz
    assert lcs4[0,6] == 0 # xyyy
    assert lcs4[0,7] == 0 # xyyz
    assert lcs4[0,8] == 0 # xyzz
    assert lcs4[0,9] == 0 # xzzz
    assert lcs4[0,10] == 3/S(8) # yyyy
    assert lcs4[0,11] == 0 # yyyz
    assert lcs4[0,12] == -3*sqrt(3)/sqrt(35) # yyzz
    assert lcs4[0,13] == 0 # yzzz
    assert lcs4[0,14] == 1 # zzzz

def test_lcs_h():
    lcs4 = get_poly_conversion(5)
    assert lcs4[0,0] == 0 # xxxxx
    assert lcs4[0,1] == 0 # xxxxy
    assert lcs4[0,2] == 5/S(8) # xxxxz
    assert lcs4[0,3] == 0 # xxxyy
    assert lcs4[0,4] == 0 # xxxyz
    assert lcs4[0,5] == 0 # xxxzz
    assert lcs4[0,6] == 0 # xxyyy
    assert lcs4[0,7] == sqrt(15)/S(4)/sqrt(7) # xxyyz
    assert lcs4[0,8] == 0 # xxyzz
    assert lcs4[0,9] == -5/sqrt(21) # xxzzz
    assert lcs4[0,10] == 0 # xyyyy
    assert lcs4[0,11] == 0 # xyyyz
    assert lcs4[0,12] == 0 # xyyzz
    assert lcs4[0,13] == 0 # xyzzz
    assert lcs4[0,14] == 0 # xzzzz
    assert lcs4[0,15] == 0 # yyyyy
    assert lcs4[0,16] == 5/S(8) # yyyyz
    assert lcs4[0,17] == 0 # yyyzz
    assert lcs4[0,18] == -5/sqrt(21) # yyzzz
    assert lcs4[0,19] == 0 # yzzzz
    assert lcs4[0,20] == 1 # zzzzz
