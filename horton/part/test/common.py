# -*- coding: utf-8 -*-
# Horton is a development platform for electronic structure methods.
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
from glob import glob

from horton import *


__all__ = [
    'get_proatomdb_cp2k', 'get_proatomdb_hf_sto3g',
    'get_proatomdb_hf_lan', 'get_fake_co', 'get_fake_pseudo_oo',
    'check_names', 'check_proatom_splines',
]


def get_proatomdb_cp2k():
    '''Return a proatomdb of pseudo oxygens and one silicon for testing purposes'''
    fns = glob(context.get_fn('test/atom_*.cp2k.out'))
    return ProAtomDB.from_files(fns)


def get_proatomdb_hf_sto3g():
    '''Return a proatomdb of H and O at hf/sto-3g for testing purposes'''
    fns = glob(context.get_fn('test/atom_???_???_hf_sto3g.fchk'))
    return ProAtomDB.from_files(fns)


def get_proatomdb_hf_lan():
    '''Return a proatomdb of H, O, Si at hf/LANL2MB for testing purposes'''
    fns = glob(context.get_fn('test/atom_???_???_hf_lan.fchk'))
    return ProAtomDB.from_files(fns)


def get_fake_co():
    # Define system
    numbers = np.array([6, 8])
    coordinates = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.132]])
    sys = System(coordinates, numbers)

    # Load some pro-atoms
    proatomdb = ProAtomDB.from_refatoms(numbers=[6, 8], max_kation=1, max_anion=1)
    proatomdb.compact(0.02)

    # Make fake cube data
    origin = np.array([-3.0, -3.0, -3.0])
    rvecs = np.identity(3, float)*0.2
    shape = np.array([30, 30, 30+11])
    ugrid = UniformGrid(origin, rvecs, shape, np.ones(3, int))

    moldens = np.zeros(ugrid.shape)
    setup = [
        (0, {+1: 0.5, 0: 0.4, -1: 0.1}),
        (1, {+1: 0.1, 0: 0.4, -1: 0.5}),
    ]
    for i, lico in setup:
        n = sys.numbers[i]
        c = sys.coordinates[i]
        spline = proatomdb.get_spline(n, lico)
        ugrid.eval_spline(spline, c, moldens)

    return sys, ugrid, moldens, proatomdb


def get_fake_pseudo_oo():
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
    ugrid = UniformGrid(origin, rvecs, shape, np.ones(3, int))

    moldens = np.zeros(ugrid.shape)
    setup = [
        (0, {+1: 0.5, 0: 0.4, -1: 0.1}),
        (1, {+1: 0.1, 0: 0.4, -1: 0.5}),
    ]
    for i, lico in setup:
        n = sys.numbers[i]
        c = sys.coordinates[i]
        spline = proatomdb.get_spline(n, lico)
        ugrid.eval_spline(spline, c, moldens)

    return sys, ugrid, moldens, proatomdb


def check_names(names, part):
    for name in names:
        assert name in part.cache


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
