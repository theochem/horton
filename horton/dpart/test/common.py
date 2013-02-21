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


__all__ = ['get_proatomdb_ref', 'get_proatomdb_cp2k', 'get_proatomdb_hf_sto3g']


def get_proatomdb_ref(numbers, max_kation, max_anion):
    '''Return a proatomdb for testing purposes'''
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rtf, int1d, 110), random_rotate=False, keep_subgrids=1)
    return ProAtomDB.from_refatoms(atgrid, numbers, max_kation, max_anion)


def get_proatomdb_cp2k():
    '''Return a proatomdb of pseudo oxygens and one silicon for testing purposes'''
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rtf, int1d, 110), random_rotate=False, keep_subgrids=1)
    fns = glob(context.get_fn('test/atom_*.cp2k.out'))
    return ProAtomDB.from_files(fns, atgrid)


def get_proatomdb_hf_sto3g():
    '''Return a proatomdb of H and O at hf/sto-3g for testing purposes'''
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 100)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rtf, int1d, 110), random_rotate=False, keep_subgrids=1)
    fns = glob(context.get_fn('test/atom_???_???_hf_sto3g.fchk'))
    return ProAtomDB.from_files(fns, atgrid)
