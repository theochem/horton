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

from horton.dpart.test.common import get_proatomdb_ref, get_proatomdb_cp2k


__all__ = ['get_fake_co', 'get_fake_pseudo_oo']


def get_fake_co():
    # Define system
    numbers = np.array([6, 8])
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.132]])
    sys = System(coordinates, numbers)

    # Load some pro-atoms
    proatomdb = get_proatomdb_ref([6, 8], max_kation=1, max_anion=1)

    # Make fake cube data
    origin = np.array([-3.0, -3.0, -3.0])
    rvecs = np.identity(3, float)*0.1
    shape = np.array([60, 60, 60+22])
    ui_grid = UniformIntGrid(origin, rvecs, shape, np.ones(3, int))

    mol_dens = np.zeros(ui_grid.shape)
    tmp = np.zeros(ui_grid.shape)
    setup = [
        (0, +1, 0.5),
        (0,  0, 0.4),
        (0, -1, 0.1),
        (1, +1, 0.1),
        (1,  0, 0.4),
        (1, -1, 0.5),
    ]
    for i, charge, frac in setup:
        n = sys.numbers[i]
        c = sys.coordinates[i]
        # TODO: can be made more efficient by scaling the spline function
        spline = proatomdb.get_spline(n, charge)
        tmp[:] = 0.0
        ui_grid.eval_spline(spline, c, tmp)
        tmp *= frac
        mol_dens += tmp

    return sys, ui_grid, mol_dens, proatomdb


def get_fake_pseudo_oo(smooth=False):
    # Define system
    numbers = np.array([8, 8])
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.132]])
    pseudo_numbers = np.array([6, 6])
    sys = System(coordinates, numbers, pseudo_numbers=pseudo_numbers)

    # Load some pro-atoms
    proatomdb = get_proatomdb_cp2k()

    # Make fake cube data
    origin = np.array([-3.0, -3.0, -3.0])
    rvecs = np.identity(3, float)*0.1
    shape = np.array([60, 60, 60+22])
    ui_grid = UniformIntGrid(origin, rvecs, shape, np.ones(3, int))

    mol_dens = np.zeros(ui_grid.shape)
    tmp = np.zeros(ui_grid.shape)
    setup = [
        (0, +1, 0.5),
        (0,  0, 0.4),
        (0, -1, 0.1),
        (1, +1, 0.1),
        (1,  0, 0.4),
        (1, -1, 0.5),
    ]
    for i, charge, frac in setup:
        n = sys.numbers[i]
        c = sys.coordinates[i]
        # TODO: can be made more efficient by scaling the spline function
        spline = proatomdb.get_spline(n, charge)
        tmp[:] = 0.0
        ui_grid.eval_spline(spline, c, tmp)
        tmp *= frac
        mol_dens += tmp

    return sys, ui_grid, mol_dens, proatomdb
