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


def check_volume_elements(rtf):
    eps = 1e-4
    npoint = 10
    r1 = rtf.get_radii(npoint, -0.5*eps)
    r2 = rtf.get_radii(npoint, +0.5*eps)
    g_ana = rtf.get_volume_elements(npoint)
    g_numer = (r2 - r1)/eps
    assert abs(g_ana - g_numer).max() < 1e-5


def test_volume_elements_log():
    rtf = LogRTransform(None, 0.1, 1e1, 100)
    assert rtf.rmin == 0.1
    assert rtf.rmax == 1e1
    assert rtf.npoint == 100
    assert rtf.alpha > 0
    check_volume_elements(rtf)
