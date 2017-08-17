# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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


import os
import shutil
import tempfile
import numpy as np

from glob import glob
from contextlib import contextmanager

from .. proatomdb import ProAtomDB, ProAtomRecord
from horton.grid import LinearRTransform, PowerRTransform, ExpRTransform, RadialGrid


__all__ = [
    'get_fn', 'load_molecule_npz', 'load_atoms_npz', 'check_names', 'check_proatom_splines',
]


def get_fn(fn):
    cur_pth = os.path.split(__file__)[0]
    return "{0}/cached/{1}".format(cur_pth, fn)


def load_molecule_npz(filename, spin_dens=False):
    # get file path
    filepath = get_fn(filename)
    # load npz file
    with np.load(filepath, mmap_mode=None) as npz:
        dens = npz['dens']
        points = npz['points']
        numbers = npz['numbers']
        coordinates = npz['coordinates']
        pseudo_numbers = npz['pseudo_numbers']
        if spin_dens:
            spindens = npz['spin_dens']
            return coordinates, numbers, pseudo_numbers, dens, spindens, points
    return coordinates, numbers, pseudo_numbers, dens, points


def get_atoms_npz(numbers, max_cation, max_anion, rtf_type, level):
    filepaths = []
    for number in numbers:
        nelectrons = number - np.arange(max_anion, max_cation + 1)
        for nelec in nelectrons:
            if level is None:
                filename = get_fn('atom_Z%.2i_N%.2i_%s.npz' % (number, nelec, rtf_type))
            else:
                filename = get_fn('atom_%s_Z%.2i_N%.2i_%s.npz' % (level, number, nelec, rtf_type))
            if os.path.isfile(filename):
                filepaths.append(filename)
    return filepaths


def load_atoms_npz(numbers, max_cation, max_anion, rtf_type='pow', level=None):
    # radial transformation available
    rtf_classes = {'lin': LinearRTransform, 'exp': ExpRTransform, 'pow': PowerRTransform}
    # get filepath of atoms npz
    filepaths = get_atoms_npz(numbers, max_cation, max_anion, rtf_type, level)
    # load each file into a record
    records = []
    rtfclass = rtf_classes[rtf_type]
    for filepath in filepaths:
        with np.load(filepath, mmap_mode=None) as npz:
            number, charge, energy = int(npz['number']), int(npz['charge']), float(npz['energy'])
            dens, deriv = npz['dens'], npz['deriv']
            rgrid = RadialGrid(rtfclass(*npz['rgrid']))
            if 'pseudo_number' in npz.keys():
                pseudo_number = npz['pseudo_number']
            else:
                pseudo_number = None
            record = ProAtomRecord(number, charge, energy, rgrid, dens, deriv, pseudo_number=pseudo_number)
            records.append(record)
    return records


@contextmanager
def tmpdir(name):
    dn = tempfile.mkdtemp(name)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)


def check_names(names, part):
    for name in names:
        assert name in part.cache


def check_proatom_splines(part):
    for index in range(part.natom):
        spline = part.get_proatom_spline(index)
        grid = part.get_grid(index)
        array1 = grid.zeros()
        part.eval_spline(index, spline, array1, grid)
        array2 = grid.zeros()
        part.eval_proatom(index, array2, grid)
        assert abs(array1).max() != 0.0
        assert abs(array1 - array2).max() < 1e-5
