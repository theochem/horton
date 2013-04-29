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


import tempfile, shutil, os, h5py as h5

from horton import *
from horton.test.common import check_script
from horton.scripts.test.common import copy_files, check_files


def write_atomdb_refatoms(tmpdir):
    int1d = SimpsonIntegrator1D()
    rtf = ExpRTransform(1e-3, 1e1, 30)
    atgrid = AtomicGrid(0, np.zeros(3, float), (rtf, int1d, 26))
    padb = ProAtomDB.from_refatoms(atgrid, numbers=[8,14], max_kation=3, max_anion=3)
    padb.to_file(os.path.join(tmpdir, 'atoms.h5'))


def check_script_jbw_coarse(scheme):
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_cpart.test_script_jbw_coarse_%s' % scheme)
    try:
        fn_cube = 'jbw_coarse_aedens.cube'
        copy_files(tmpdir, [fn_cube])
        write_atomdb_refatoms(tmpdir)
        check_script('horton-cpart.py %s %s atoms.h5' % (fn_cube, scheme), tmpdir)
        fn_h5 = '%s.h5' % fn_cube
        check_files(tmpdir, [fn_h5])
        with h5.File(os.path.join(tmpdir, fn_h5)) as f:
            assert 'cpart' in f
            assert scheme + '_r1' in f['cpart']
    finally:
        shutil.rmtree(tmpdir)


def test_script_jbw_coarse_h():
    # other schemes don't work because the cube file is too crappy.
    check_script_jbw_coarse('h')
