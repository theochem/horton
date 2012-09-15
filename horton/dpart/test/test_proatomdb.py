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


import numpy as np, h5py as h5
import tempfile, os
from horton import *


def test_from_scratch_simple():
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 1e1, 100)
    atgrid = AtomicGrid(np.zeros(3, float), rtf, 110, 100, keep_subgrids=1)
    pro_atomdb = ProAtomDB.from_scratch([HartreeFock()], '3-21G', atgrid, [1,6], qmin=-2, qmax=3)
    keys = sorted(pro_atomdb._records.keys())
    assert keys == [(1, 1), (1, 2), (1, 3), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7), (6, 8)]


def test_io_group():
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 1e1, 100)
    atgrid = AtomicGrid(np.zeros(3, float), rtf, 110, 100, keep_subgrids=1)
    pro_atomdb = ProAtomDB.from_scratch([HartreeFock()], '3-21G', atgrid, [1,6], qmin=-1, qmax=1)
    keys = sorted(pro_atomdb._records.keys())
    assert keys == [(1, 1), (1, 2), (6, 5), (6, 6), (6, 7)]

    with h5.File('horton.dpart.test.test_proatomdb.test_io_group', driver='core', backing_store=False) as f:
        pro_atomdb.to_file(f)
        bis = ProAtomDB.from_file(f)
        keys = sorted(bis._records.keys())
        assert keys == [(1, 1), (1, 2), (6, 5), (6, 6), (6, 7)]
        avr1 = pro_atomdb._records[(1,2)]
        avr2 = bis._records[(1,2)]
        assert (avr1 == avr2).all()

def test_io_filename():
    rtf = LogRTransform(TrapezoidIntegrator1D(), 1e-3, 1e1, 100)
    atgrid = AtomicGrid(np.zeros(3, float), rtf, 110, 100, keep_subgrids=1)
    pro_atomdb = ProAtomDB.from_scratch([HartreeFock()], '3-21G', atgrid, [1,6], qmin=0, qmax=1)
    keys = sorted(pro_atomdb._records.keys())
    assert keys == [(1, 1), (6, 5), (6, 6)]

    tmpdir = tempfile.mkdtemp('horton.dpart.test.test_proatomdb.test_io_filename')
    filename = '%s/test.h5' % tmpdir
    try:
        pro_atomdb.to_file(filename)
        bis = ProAtomDB.from_file(filename)
        keys = sorted(bis._records.keys())
        assert keys == [(1, 1), (6, 5), (6, 6)]
        avr1 = pro_atomdb._records[(6,5)]
        avr2 = bis._records[(6,5)]
        assert (avr1 == avr2).all()
    finally:
        if os.path.isfile(filename):
            os.remove(filename)
        os.rmdir(tmpdir)
