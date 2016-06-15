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


import numpy as np, h5py as h5
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import

from horton.test.common import get_random_cell, compare_symmetries


def test_symmetry_attrs():
    generators = np.random.uniform(-1, 1, (5, 3, 4))
    fracs = np.random.uniform(0, 1, (4, 3))
    numbers = np.array([1, 6, 6, 1])
    cell = get_random_cell(10.0, 3)
    s = Symmetry('boo', generators, fracs, numbers, cell)
    assert (s.labels == ['H0', 'C1', 'C2', 'H3']).all()
    s = Symmetry('boo', generators, fracs, numbers, cell, ['q', 'w', 'e', 'r'])
    assert s.name == 'boo'
    assert (s.generators == generators).all()
    assert s.natom == 4
    assert s.fracs is s.fracs
    assert s.numbers is numbers
    assert s.cell is cell


def get_fake_symmetry():
    g1 = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
    ])
    g2 = np.array([
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 1, 0],
    ], dtype=float)
    generators = [g1, g2]
    fracs = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 0.0, 0.5]])
    numbers = np.array([6, 1, 8])
    cell = Cell(np.diag([10.0, 12.0, 14.0]))
    return Symmetry('fried eggs', generators, fracs, numbers, cell)

def test_generate():
    s = get_fake_symmetry()
    coordinates, numbers, links = s.generate()
    assert coordinates.shape == (4, 3)
    assert numbers.shape == (4,)
    assert links.shape == (4, 2)
    assert abs(coordinates[0] - [0.0, 0.0, 0.0]).max() < 1e-10
    assert abs(coordinates[1] - [5.0, 0.0, 0.0]).max() < 1e-10
    assert abs(coordinates[2] - [0.0, 6.0, 0.0]).max() < 1e-10
    assert abs(coordinates[3] - [0.0, 0.0, 7.0]).max() < 1e-10
    assert (numbers == [6, 1, 1, 8]).all()
    assert (links == [[0, 0], [1, 0], [1, 1], [2, 0]]).all()


def test_identify():
    mol1 = IOData.from_file(context.get_fn('test/lta_gulp.cif'))
    mol2 = IOData.from_file(context.get_fn('test/aelta.cube'))
    sym = mol1.symmetry
    links = sym.identify(mol2.coordinates, mol1.cell) # The cell parameters in aelta.cube are rubbish.
    for i in xrange(mol2.natom):
        if mol2.numbers[i] == 14:
            assert links[i,0] == 0


def test_hdf5():
    s0 = get_fake_symmetry()
    with h5.File('horton.test.test_symmetry.test_hdf5.h5', driver='core', backing_store=False) as f:
        s0.to_hdf5(f)
        s1 = Symmetry.from_hdf5(f)
        compare_symmetries(s0, s1)


def test_symmetry_error():
    mol1 = IOData.from_file(context.get_fn('test/lta_gulp.cif'))
    mol2 = IOData.from_file(context.get_fn('test/lta_iza.cif'))
    with assert_raises(SymmetryError):
        mol1.symmetry.identify(mol2.coordinates, mol2.cell)
