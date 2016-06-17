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


from glob import glob

from horton import ProAtomDB, context


__all__ = [
    'get_proatomdb_cp2k', 'get_proatomdb_hf_sto3g',
    'get_proatomdb_hf_lan', 'check_names', 'check_proatom_splines',
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


def check_names(names, part):
    for name in names:
        assert name in part.cache


def check_proatom_splines(part):
    for index in xrange(part.natom):
        spline = part.get_proatom_spline(index)
        grid = part.get_grid(index)
        array1 = grid.zeros()
        part.eval_spline(index, spline, array1, grid)
        array2 = grid.zeros()
        part.eval_proatom(index, array2, grid)
        assert abs(array1).max() != 0.0
        assert abs(array1 - array2).max() < 1e-5
