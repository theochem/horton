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

from horton.grid.uniform import UniformIntGrid
from horton.espfit.cext import setup_esp_cost_cube


__all__ = ['ESPCost']


class ESPCost(object):
    def __init__(self, system, grid, vref, weights, rcut=20, alpha_scale=3.0, gcut_scale=1.1):
        alpha = alpha_scale / rcut
        gcut = gcut_scale * alpha
        self._A = np.zeros((system.natom, system.natom), float)
        self._B = np.zeros(system.natom, float)
        self._C = np.zeros((), float)
        if isinstance(grid, UniformIntGrid):
            setup_esp_cost_cube(grid.origin, grid.grid_cell, grid.shape,
                grid.get_cell(), vref, weights, system.coordinates, self._A, self._B, self._C, rcut, alpha, gcut)
        else:
            raise NotImplementedError

    def value(self, charges):
        return 0.5*np.dot(charges, np.dot(self._A, charges)) + np.dot(charges, self._B) + 0.5*self._C

    def gradient(self, charges):
        return np.dot(self._A, charges) + self._B
