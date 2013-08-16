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


import h5py as h5, argparse, numpy as np, os
from nose.tools import assert_raises

from horton import *
from horton.scripts.common import *
from horton.test.common import tmpdir


def test_get_output_filename():
    assert get_output_filename('some.file', 'boo') == 'some_boo.h5'
    assert get_output_filename('some_file', 'boo') == 'some_file_boo.h5'
    assert get_output_filename('some.file', 'boo', 'bar') == 'bar'
    assert get_output_filename('some_file', 'boo', 'bar.h5') == 'bar.h5'


def test_iter_elements():
    assert list(iter_elements('1,2')) == [1, 2]
    assert list(iter_elements('H,2')) == [1, 2]
    assert list(iter_elements('H,He')) == [1, 2]
    assert list(iter_elements('1,He')) == [1, 2]
    assert list(iter_elements('1-6')) == [1, 2, 3, 4, 5, 6]
    assert list(iter_elements('1-C')) == [1, 2, 3, 4, 5, 6]
    assert list(iter_elements('H-C')) == [1, 2, 3, 4, 5, 6]
    assert list(iter_elements('H-6')) == [1, 2, 3, 4, 5, 6]
    assert list(iter_elements('6-8')) == [6, 7, 8]
    assert list(iter_elements('2,6-8')) == [2,6, 7, 8]
    assert list(iter_elements('2,6-8,10')) == [2, 6, 7, 8, 10]
    assert list(iter_elements('10,6-8,2')) == [10 ,6, 7, 8, 2]
    assert list(iter_elements('6-8,2')) == [6, 7, 8, 2]
    with assert_raises(ValueError):
        list(iter_elements('8-6,2')) == [2]


def test_reduce_data0():
    from horton import UniformGrid
    data = np.random.normal(0, 1, (10, 20, 30))
    grid_rvecs = np.identity(3, float)*0.1
    ugrid = UniformGrid(np.array([0.5, 0.3, -0.1]), grid_rvecs, np.array(data.shape), np.array([0, 1, 0]))

    data1, ugrid1 = reduce_data(data, ugrid, 10, 0)
    assert data1.shape == (1, 2, 3)
    assert (data1 == data[::10,::10,::10]).all()
    assert (ugrid1.origin == ugrid.origin).all()
    assert abs(ugrid1.grid_cell.rvecs - ugrid.grid_cell.rvecs*10).max() < 1e-10
    assert (ugrid1.shape == [1, 2, 3]).all()
    assert (ugrid1.pbc == ugrid.pbc).all()

    data2, ugrid2 = reduce_data(data, ugrid, 5, 0)
    assert data2.shape == (2, 4, 6)
    assert (data2 == data[::5,::5,::5]).all()
    assert (ugrid2.origin == ugrid.origin).all()
    assert abs(ugrid2.grid_cell.rvecs - ugrid.grid_cell.rvecs*5).max() < 1e-10
    assert (ugrid2.shape == [2, 4, 6]).all()
    assert (ugrid2.pbc == ugrid.pbc).all()


def test_reduce_data1():
    from horton import UniformGrid
    data = np.random.normal(0, 1, (11, 21, 31))
    grid_rvecs = np.identity(3, float)*0.1
    ugrid = UniformGrid(np.array([0.3, 0.2, -0.1]), grid_rvecs, np.array(data.shape), np.array([1, 1, 0]))

    data1, ugrid1 = reduce_data(data, ugrid, 10, 1)
    assert data1.shape == (1, 2, 3)
    assert (data1 == data[:-1:10,:-1:10,:-1:10]).all()
    assert (ugrid1.origin == ugrid.origin).all()
    assert abs(ugrid1.grid_cell.rvecs - ugrid.grid_cell.rvecs*10).max() < 1e-10
    assert (ugrid1.shape == [1, 2, 3]).all()
    assert (ugrid1.pbc == ugrid.pbc).all()

    data2, ugrid2 = reduce_data(data, ugrid, 5, 1)
    assert data2.shape == (2, 4, 6)
    assert (data2 == data[:-1:5,:-1:5,:-1:5]).all()
    assert (ugrid2.origin == ugrid.origin).all()
    assert abs(ugrid2.grid_cell.rvecs - ugrid.grid_cell.rvecs*5).max() < 1e-10
    assert (ugrid2.shape == [2, 4, 6]).all()
    assert (ugrid2.pbc == ugrid.pbc).all()


def test_parse_pbc():
    assert (parse_pbc('111') == [1, 1, 1]).all()
    assert (parse_pbc('000') == [0, 0, 0]).all()
    assert (parse_pbc('101') == [1, 0, 1]).all()
    assert (parse_pbc('001') == [0, 0, 1]).all()


def test_parse_ugrid_1():
    rvecs = np.diag([3.0, 2.0, 1.0])*angstrom
    cell = Cell(rvecs)
    ugrid = parse_ugrid('0.1', cell)

    assert (ugrid.origin == [0, 0, 0]).all()
    assert abs(ugrid.grid_cell.rvecs - np.identity(3, float)*0.1*angstrom).max() < 1e-10
    assert (ugrid.shape == [30, 20, 10]).all()
    assert (ugrid.pbc == 1).all()


def test_parse_ugrid_2():
    with tmpdir('horton.scripts.test.test_common.test_parse_ugrid_2') as dn:
        fn_h5 = os.path.join(dn, 'test.h5')
        origin = np.random.uniform(0, 1, 3)
        grid_rvecs = np.random.uniform(0, 1, (3, 3))
        shape = np.random.randint(10, 20, 3)
        pbc = np.random.randint(0, 2, 3)
        ugrid1 = UniformGrid(origin, grid_rvecs, shape, pbc)

        with h5.File(fn_h5) as f:
            ugrid1.to_hdf5(f)

        ugrid2 = parse_ugrid('%s:/' % fn_h5, None)

        assert (ugrid2.origin == origin).all()
        assert (ugrid2.grid_cell.rvecs == grid_rvecs).all()
        assert (ugrid2.shape == shape).all()
        assert (ugrid2.pbc == pbc).all()


def test_store_args():
    with h5.File('horton.scripts.test.test_common.test_store_args.h5', driver='core', backing_store=False) as f:
        namespace = argparse.Namespace()
        namespace.a = 1
        namespace.b = 'a'
        namespace.c = None
        namespace.d = 2.5
        store_args(namespace, f)
        assert f.attrs['arg_a'] == 1
        assert f.attrs['arg_b'] == 'a'
        assert 'arg_c' not in f.attrs
        assert f.attrs['arg_d'] == 2.5
        assert 'pwd' in f.attrs
        assert 'cmdline' in f.attrs
        assert 'datetime' in f.attrs


def test_safe_open1():
    # just a silly test
    with safe_open_h5('horton.scripts.test.test_common.test_safe_open1.h5', driver='core', backing_store=False) as f:
        pass


def test_safe_open2():
    # test error handling in with clause
    with assert_raises(ValueError):
        with safe_open_h5('horton.scripts.test.test_common.test_safe_open2.h5', driver='core', backing_store=False) as f:
            raise ValueError


def test_safe_open3():
    # test error handling in file opening
    with assert_raises(ValueError):
        with safe_open_h5('horton.scripts.test.test_common.test_safe_open3.h5', driver='fubar', wait=0.1, count=3) as f:
            raise ValueError
