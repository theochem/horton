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


import csv
import h5py as h5, numpy as np

from horton.test.common import check_script, tmpdir
from horton.scripts.test.common import check_files
from horton.scripts.hdf2csv import iter_datasets


def fill_hdf5(f):
    f['test1'] = np.random.uniform(-1, 1, 5)
    f['test2'] = np.random.uniform(-1, 1, (5, 2))
    f['test0'] = np.random.uniform(-1, 1)
    g = f.create_group('bar')
    g['foo'] = np.random.uniform(-1, 1, (5, 2))
    f['zzz'] = np.random.uniform(-1, 1, (5, 2, 5))


def test_iter_datasets():
    with h5.File('horton.scripts.test.test_hdf2csv.test_iter_datasets', driver='core', backing_store=False) as f:
        fill_hdf5(f)
        l = list(iter_datasets(f))
        assert len(l) == 5
        assert l[0][0] == 'bar/foo'
        assert l[0][1] == f['bar/foo']
        assert l[1][0] == 'test0'
        assert l[1][1] == f['test0']
        assert l[2][0] == 'test1'
        assert l[2][1] == f['test1']
        assert l[3][0] == 'test2'
        assert l[3][1] == f['test2']
        assert l[4][0] == 'zzz'
        assert l[4][1] == f['zzz']


def test_script():
    with tmpdir('horton.scripts.test.test_hdf2csv.test_script') as dn:
        with h5.File('%s/test.h5' % dn) as f:
            fill_hdf5(f)
        check_script('horton-hdf2csv.py test.h5:/ test.csv', dn)
        check_files(dn, ['test.h5', 'test.csv'])
        with open('%s/test.csv' % dn) as f:
            r = csv.reader(f)
            rows = [row for row in r]
        with h5.File('%s/test.h5' % dn) as f:
            assert rows[0][0] == 'Converted data from test.h5:/'
            assert len(rows[1]) == 0
            assert len(rows[2]) == 2
            assert rows[2][0] == 'Dataset'
            assert rows[2][1] == 'bar/foo'
            assert len(rows[3]) == 3
            assert rows[3][0] == 'Shape'
            assert rows[3][1] == '5'
            assert rows[3][2] == '2'
            assert float(rows[4][0]) == f['bar/foo'][0,0]
            assert float(rows[8][1]) == f['bar/foo'][4,1]
            assert len(rows[9]) == 0
            assert len(rows[22]) == 2
            assert rows[22][0] == 'Dataset'
            assert rows[22][1] == 'test2'
            assert len(rows[23]) == 3
            assert rows[23][0] == 'Shape'
            assert rows[23][1] == '5'
            assert rows[23][2] == '2'
            assert float(rows[24][0]) == f['test2'][0,0]
            assert float(rows[28][1]) == f['test2'][4,1]
            assert rows[29] == []
            assert len(rows[30]) == 2
            assert rows[30][0] == 'Dataset'
            assert rows[30][1] == 'zzz'
            assert len(rows[31]) == 4
            assert rows[31][0] == 'Shape'
            assert rows[31][1] == '5'
            assert rows[31][2] == '2'
            assert rows[31][3] == '5'
            assert float(rows[32][0]) == f['zzz'][0,0,0]
            assert float(rows[32][1]) == f['zzz'][0,1,0]
            assert rows[32][2] == ''
            assert float(rows[32][3]) == f['zzz'][0,0,1]
            assert float(rows[32][4]) == f['zzz'][0,1,1]
            assert rows[32][5] == ''
            assert len(rows[32]) == 3*5-1
            assert float(rows[36][1]) == f['zzz'][4,1,0]
            assert float(rows[36][3*5-2]) == f['zzz'][4,1,4]
            assert rows[-1] == []
