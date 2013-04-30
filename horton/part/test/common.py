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
from glob import glob

from horton import *


__all__ = [
    'get_proatomdb_ref', 'get_proatomdb_cp2k', 'get_proatomdb_hf_sto3g',
    'get_proatomdb_hf_lan', 'get_fake_co', 'get_fake_pseudo_oo',
    'check_names', 'check_proatom_splines',
]


def get_proatomdb_ref(numbers, max_kation, max_anion):
    '''Return a proatomdb for testing purposes'''
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialIntGrid(rtf)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rgrid, 110), random_rotate=False)
    return ProAtomDB.from_refatoms(atgrid, numbers, max_kation, max_anion)


def get_proatomdb_cp2k():
    '''Return a proatomdb of pseudo oxygens and one silicon for testing purposes'''
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialIntGrid(rtf)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rgrid, 110), random_rotate=False)
    fns = glob(context.get_fn('test/atom_*.cp2k.out'))
    return ProAtomDB.from_files(fns, atgrid)


def get_proatomdb_hf_sto3g():
    '''Return a proatomdb of H and O at hf/sto-3g for testing purposes'''
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialIntGrid(rtf)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rgrid, 110), random_rotate=False)
    fns = glob(context.get_fn('test/atom_???_???_hf_sto3g.fchk'))
    return ProAtomDB.from_files(fns, atgrid)


def get_proatomdb_hf_lan():
    '''Return a proatomdb of H, O, Si at hf/LANL2MB for testing purposes'''
    rtf = ExpRTransform(5e-4, 2e1, 120)
    rgrid = RadialIntGrid(rtf)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rgrid, 110), random_rotate=False)
    fns = glob(context.get_fn('test/atom_???_???_hf_lan.fchk'))
    return ProAtomDB.from_files(fns, atgrid)


def get_fake_co():
    # Define system
    numbers = np.array([6, 8])
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.132]])
    sys = System(coordinates, numbers)

    # Load some pro-atoms
    proatomdb = get_proatomdb_ref([6, 8], max_kation=1, max_anion=1)
    proatomdb.compact(0.02)

    # Make fake cube data
    origin = np.array([-3.0, -3.0, -3.0])
    rvecs = np.identity(3, float)*0.2
    shape = np.array([30, 30, 30+11])
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
    proatomdb.compact(0.02)

    # Make fake cube data
    origin = np.array([-3.0, -3.0, -3.0])
    rvecs = np.identity(3, float)*0.2
    shape = np.array([30, 30, 30+11])
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


def check_names(names, part):
    for name in names:
        assert part.cache.has(name)


def check_proatom_splines(part):
    for index in xrange(part.system.natom):
        spline = part.get_proatom_spline(index)
        grid = part.get_grid(index)
        array1 = grid.zeros()
        part.eval_spline(index, spline, array1)
        array2 = grid.zeros()
        part.eval_proatom(index, array2)
        assert abs(array1).max() != 0.0
        assert abs(array1 - array2).max() < 1e-5
