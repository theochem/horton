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

from horton.grid.cext import LinearRTransform, CubicSpline


__all__ = ['get_cosine_spline', 'get_exp_spline']


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
