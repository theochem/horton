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

from horton import System, context
from horton.scripts.wpart import *
from horton.test.common import check_script
from horton.scripts.test.common import copy_files, check_files
from horton.part.test.common import get_proatomdb_ref, get_proatomdb_hf_sto3g


def test_parse_grid_1():
    for sgrid in 'coarse', 'medium', 'fine', 'veryfine':
        atspecs = parse_grid(sgrid, None, None)
        fn = context.get_fn('grids/%s.txt' % atspecs)
        assert os.path.isfile(fn)


def test_parse_grid_2():
    sys = System.from_file(context.get_fn('test/water_number.xyz'))
    padb = get_proatomdb_ref([1, 8], 1, 1)
    for sgrid in '6', '38', '110':
        nll = int(sgrid)
        atspecs = parse_grid(sgrid, sys, padb)
        for i in xrange(sys.natom):
            assert len(atspecs[i]) == 3
            assert atspecs[i][0] == padb.get_rgrid(sys.numbers[i]).rtransform
            assert atspecs[i][2] == nll


def write_atomdb_sto3g(tmpdir):
    padb = get_proatomdb_hf_sto3g()
    padb.to_file(os.path.join(tmpdir, 'atoms.h5'))


def check_script_water_sto3g(scheme):
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_wpart.test_script_water_sto3g')
    try:
        fn_fchk = 'water_sto3g_hf_g03.fchk'
        copy_files(tmpdir, [fn_fchk])
        write_atomdb_sto3g(tmpdir)
        check_script('horton-wpart.py water_sto3g_hf_g03.fchk %s atoms.h5' % scheme, tmpdir)
        fn_h5 = 'water_sto3g_hf_g03.fchk.h5'
        check_files(tmpdir, [fn_h5])
        with h5.File(os.path.join(tmpdir, fn_h5)) as f:
            assert 'wpart' in f
            assert scheme in f['wpart']
    finally:
        shutil.rmtree(tmpdir)


def test_script_water_sto3g_h():
    check_script_water_sto3g('h')


def test_script_water_sto3g_hi():
    check_script_water_sto3g('hi')


def test_script_water_sto3g_he():
    check_script_water_sto3g('he')
