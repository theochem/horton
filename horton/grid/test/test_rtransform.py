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

def test_string_log_rtransform():
    rtf1 = LogRTransform(None, np.random.uniform(1e-5, 5e-5), np.random.uniform(1, 5), 111)
    s = rtf1.to_string()
    rtf2 = BaseRTransform.from_string(s, None)
    assert rtf1.rmin == rtf2.rmin
    assert rtf1.rmax == rtf2.rmax
    assert rtf1.npoint == rtf2.npoint
    assert rtf1.alpha == rtf2.alpha

    try:
        rtf3 = BaseRTransform.from_string('Fubar A 5', None)
        assert False
    except TypeError:
        pass

    try:
        rtf3 = BaseRTransform.from_string('LogRTransform A 5', None)
        assert False
    except ValueError:
        pass

    try:
        rtf3 = BaseRTransform.from_string('LogRTransform A 5 .1', None)
        assert False
    except ValueError:
        pass

    rtf3 = BaseRTransform.from_string('LogRTransform 1.0 12.15643216847 5', None)
    assert rtf3.rmin == 1.0
    assert rtf3.rmax == 12.15643216847
    assert rtf3.npoint == 5
