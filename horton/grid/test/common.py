# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
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


import numpy as np
from horton import *


__all__ = ['get_cosine_spline', 'get_exp_spline', 'get_random_cell']


def get_cosine_spline():
    # Construct a simple spline for the function cos(x)+1 in the range 0,pi
    rtf = LinearRTransform(0.0, np.pi, 100)
    x = rtf.get_radii()
    y = np.cos(x)+1
    d = -np.sin(x)
    return CubicSpline(y, d, rtf)


def get_exp_spline():
    rtf = LinearRTransform(0.0, 20.0, 100)
    x = rtf.get_radii()
    y = np.exp(-0.2*x)
    d = -0.2*np.exp(-0.2*x)
    return CubicSpline(y, d, rtf)


def get_random_cell(a, nvec):
    if nvec == 0:
        return Cell(None)
    while True:
        rvecs = np.random.uniform(0, a, (nvec,3))
        cell = Cell(rvecs)
        if cell.volume > a**3*0.1:
            return cell
