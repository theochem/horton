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


import tempfile, shutil, os, csv
import h5py as h5, numpy as np

from horton.test.common import check_script
from horton.scripts.test.common import check_files
from horton.scripts.hdf2csv import iter_datasets


def fill_hdf5(f):
    f['test1'] = np.random.uniform(-1, 1, 5)
    f['test2'] = np.random.uniform(-1, 1, (5, 2))
    f['test0'] = np.random.uniform(-1, 1)
    g = f.create_group('bar')
    g['foo'] = np.random.uniform(-1, 1, (5, 2))


def test_iter_datasets():
    with h5.File('horton.scripts.test.test_hdf2csv.test_iter_datasets', driver='core', backing_store=False) as f:
        fill_hdf5(f)
        l = list(iter_datasets(f))
        assert len(l) == 4
        assert l[0][0] == 'bar/foo'
        assert l[0][1] == f['bar/foo']
        assert l[1][0] == 'test0'
        assert l[1][1] == f['test0']
        assert l[2][0] == 'test1'
        assert l[2][1] == f['test1']
        assert l[3][0] == 'test2'
        assert l[3][1] == f['test2']


def test_script():
    tmpdir = tempfile.mkdtemp('horton.scripts.test.test_hdf2csv.test_script')
    try:
        with h5.File('%s/test.h5' % tmpdir) as f:
            fill_hdf5(f)
        check_script('horton-hdf2csv.py test.h5:/ test.csv', tmpdir)
        check_files(tmpdir, ['test.h5', 'test.csv'])
        os.system('cat %s/test.csv' % tmpdir)
        with open('%s/test.csv' % tmpdir) as f:
            r = csv.reader(f)
            rows = [row for row in r]
        with h5.File('%s/test.h5' % tmpdir) as f:
            assert rows[0][0] == 'Converted data from test.h5:/'
            assert len(rows[1]) == 0
            assert rows[2][0] == 'Dataset'
            assert rows[2][1] == 'bar/foo'
            assert float(rows[3][0]) == f['bar/foo'][0,0]
            assert float(rows[7][1]) == f['bar/foo'][4,1]
            assert len(rows[8]) == 0
            assert rows[19][0] == 'Dataset'
            assert rows[19][1] == 'test2'
            assert float(rows[20][0]) == f['test2'][0,0]
            assert float(rows[24][1]) == f['test2'][4,1]
            assert rows[-1] == []
    finally:
        shutil.rmtree(tmpdir)
