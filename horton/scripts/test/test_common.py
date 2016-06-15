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


import h5py as h5, argparse, numpy as np, os
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import
from horton.scripts.common import *
from horton.test.common import tmpdir


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
    assert abs(ugrid1.grid_rvecs - ugrid.grid_rvecs*10).max() < 1e-10
    assert (ugrid1.shape == [1, 2, 3]).all()
    assert (ugrid1.pbc == ugrid.pbc).all()

    data2, ugrid2 = reduce_data(data, ugrid, 5, 0)
    assert data2.shape == (2, 4, 6)
    assert (data2 == data[::5,::5,::5]).all()
    assert (ugrid2.origin == ugrid.origin).all()
    assert abs(ugrid2.grid_rvecs - ugrid.grid_rvecs*5).max() < 1e-10
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
    assert abs(ugrid1.grid_rvecs - ugrid.grid_rvecs*10).max() < 1e-10
    assert (ugrid1.shape == [1, 2, 3]).all()
    assert (ugrid1.pbc == ugrid.pbc).all()

    data2, ugrid2 = reduce_data(data, ugrid, 5, 1)
    assert data2.shape == (2, 4, 6)
    assert (data2 == data[:-1:5,:-1:5,:-1:5]).all()
    assert (ugrid2.origin == ugrid.origin).all()
    assert abs(ugrid2.grid_rvecs - ugrid.grid_rvecs*5).max() < 1e-10
    assert (ugrid2.shape == [2, 4, 6]).all()
    assert (ugrid2.pbc == ugrid.pbc).all()


def test_parse_pbc():
    assert (parse_pbc('111') == [1, 1, 1]).all()
    assert (parse_pbc('000') == [0, 0, 0]).all()
    assert (parse_pbc('101') == [1, 0, 1]).all()
    assert (parse_pbc('001') == [0, 0, 1]).all()


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


def test_check_output():
    with tmpdir('horton.scripts.test.test_common.test_check_output') as dn:
        assert not check_output('%s/foo.h5' % dn, '/', True)
        assert not check_output('%s/foo.h5' % dn, '/bar', True)
        assert not check_output('%s/foo.h5' % dn, '/', False)
        assert not check_output('%s/foo.h5' % dn, '/bar', False)
        with h5.File('%s/foo.h5' % dn) as f:
            f.create_group('bork')
        assert not check_output('%s/foo.h5' % dn, '/', True)
        assert not check_output('%s/foo.h5' % dn, '/bork', True)
        assert not check_output('%s/foo.h5' % dn, '/bar', True)
        assert check_output('%s/foo.h5' % dn, '/', False)
        assert not check_output('%s/foo.h5' % dn, '/bork', False)
        assert not check_output('%s/foo.h5' % dn, '/bar', False)
        with h5.File('%s/foo.h5' % dn) as f:
            f['bork']['a'] = np.array([1, 2, 3])
        assert not check_output('%s/foo.h5' % dn, '/', True)
        assert not check_output('%s/foo.h5' % dn, '/bork', True)
        assert not check_output('%s/foo.h5' % dn, '/bar', True)
        assert check_output('%s/foo.h5' % dn, '/', False)
        assert check_output('%s/foo.h5' % dn, '/bork', False)
        assert not check_output('%s/foo.h5' % dn, '/bar', False)


def test_write_script_output():
    class Foo:
        pass
    def test_h5(fn_h5):
        with h5.File(fn_h5, 'r') as f:
            assert sorted(f.keys()) == ['a', 'c']
            assert f['a'].keys() == ['b']
            assert (f['a/b'][:] == np.array([1, 2, 3])).all()
            assert (f['c'][:] == np.array([0.1, 0.2, 0.3])).all()
            assert len(f.attrs) == 5
    with tmpdir('horton.scripts.test.test_common.test_write_script_output') as dn:
        fn_h5 = '%s/foo.h5' % dn
        results = {'a': {'b': np.array([1, 2, 3])}, 'c': np.array([0.1, 0.2, 0.3])}
        args = Foo()
        args.bar = 'egg'
        args.bork = 'kwak'
        write_script_output(fn_h5, '/', results, args)
        test_h5(fn_h5)
        write_script_output(fn_h5, '/', results, args)
        test_h5(fn_h5)
        with h5.File(fn_h5) as f:
            f['d'] = 5
        write_script_output(fn_h5, '/', results, args)
        test_h5(fn_h5)
        with h5.File(fn_h5) as f:
            for key in f.keys():
                del f[key]
        write_script_output(fn_h5, '/', results, args)
        test_h5(fn_h5)
