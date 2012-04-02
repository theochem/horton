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

from harmonics import *

def test_lcs():
    # Comparison with some results from:
    # Schlegel, H. B.;  Frisch, M. J. Int. J. Quantum Chem. 1995, 54, 83-87.

    lcs2 = get_poly_conversion(2)
    assert lcs2[0,0] == -1/S(2) # xx
    assert lcs2[0,1] == 0 # xy
    assert lcs2[0,2] == 0 # xz
    assert lcs2[0,3] == -1/S(2) # yy
    assert lcs2[0,4] == 0 # yz
    assert lcs2[0,5] == 1 # zz

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
