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


from horton.grid.base import IntGrid
from horton.log import log


__all__ = ['RectangleGrid']


class RectangleGrid(IntGrid):
    def __init__(self, center, axis0, axis1, l0, h0, l1, h1):
        if h0 <= l0:
            raise ValueError('l0 should be lower than h0.')
        if h1 <= l1:
            raise ValueError('l1 should be lower than h1.')
        npoint = (h0-l0+1)*(h1-l1+1)
        points = np.zeros((npoint, 3), float)
        weights = np.empty(npoint, float)
        weights[:] = np.sqrt(
            (np.linalg.norm(axis0)*np.linalg.norm(axis1))**2 -
            np.dot(axis0, axis1)**2
        )
        counter = 0
        for i0 in xrange(l0, h0+1):
            for i1 in xrange(l1, h1+1):
                points[counter] = center + axis0*i0 + axis1*i1
                counter += 1
        IntGrid.__init__(self, points, weights)
