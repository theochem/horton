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


import numpy as np, h5py as h5, tempfile, os
from horton import *


class Example(JustOnceClass):
    def __init__(self):
        JustOnceClass.__init__(self)
        self.counter = 0

    @just_once
    def inc(self):
        self.counter += 1

    def inc_bis(self):
        self.counter += 1


def test_just_once():
    e = Example()
    assert e.counter == 0
    e.inc()
    assert e.counter == 1
    e.inc()
    assert e.counter == 1
    e.invalidate()
    assert e.counter == 1
    e.inc()
    assert e.counter == 2
    e.inc()
    assert e.counter == 2
    e.inc_bis()
    assert e.counter == 3


def test_cache():
    c = Cache()
    c.dump('foo', 5)
    assert c.load('foo') == 5
    c.dump('foo', 4, 6)
    assert c.load('foo', 4) == 6
    c.invalidate_all()
    assert len(c._store) == 0


def test_cache_alloc1():
    c = Cache()
    assert not c.has('bar')
    tmp, new = c.load('bar', alloc=5)
    assert new
    assert (tmp == 0).all()
    assert tmp.shape == (5,)
    assert issubclass(tmp.dtype.type, float)
    assert c.has('bar')
    tmp[3] = 1
    bis = c.load('bar')
    assert bis is tmp
    assert bis[3] == 1
    tris, new = c.load('bar', alloc=5)
    assert not new
    assert tris is tmp


def test_cache_alloc2():
    c = Cache()
    tmp, new = c.load('egg', alloc=(5,10))
    assert new
    assert (tmp == 0).all()
    assert tmp.shape == (5,10)
    assert issubclass(tmp.dtype.type, float)
    tmp[3] = 1
    bis = c.load('egg')
    assert bis is tmp
    assert (bis[3] == 1).all()
    tris, new = c.load('egg', alloc=(5,10))
    assert not new
    assert tris is tmp


def test_cache_allocation():
    c = Cache()
    assert not c.has('egg')
    tmp, new = c.load('egg', alloc=(5,10))
    assert new
    assert (tmp == 0).all()
    assert tmp.shape == (5,10)
    assert issubclass(tmp.dtype.type, float)
    assert c.has('egg')
    assert not c.has('bar')
    tmp[:] = 1.0
    c.invalidate_all()
    assert not c.has('egg')
    assert (tmp[:] == 0.0).all()
    # try to load it, while it is no longer valid
    try:
        bis = c.load('egg')
        assert False
    except KeyError:
        pass
    # properly load it anew
    bis, new = c.load('egg', alloc=(5,10))
    assert new
    assert bis is tmp # still the same array, just resetted.
    assert c.has('egg')
    # simple load should now work
    tris = c.load('egg')
    assert tris is tmp


def test_cache_default():
    c = Cache()
    assert c.load('egg', default=5)
    c.dump('egg', 5)
    assert c.load('egg') == 5
    assert c.load('egg', default=6) == 5
    c.invalidate_all()
    assert c.load('egg', default=6) == 6
    try:
        c.load('egg')
        assert False
    except KeyError:
        pass
    c.invalidate_all()
    assert c.load('egg', default=None) == None
    try:
        c.load('egg')
        assert False
    except KeyError:
        pass


def test_dense_one_body():
    from horton.matrix import DenseLinalgFactory, DenseOneBody
    lf = DenseLinalgFactory()
    c = Cache()
    op1, new = c.load('egg', alloc=(lf, 'one_body', 10))
    assert new
    assert isinstance(op1, DenseOneBody)
    assert op1.nbasis == 10
    op2 = c.load('egg')
    assert op1 is op2
    op3, new = c.load('egg', alloc=(lf, 'one_body', 10))
    assert not new
    assert op1 is op3
    # things that should not work
    try:
        op4, new = c.load('egg', alloc=(lf, 'one_body', 5))
        assert False
    except TypeError:
        pass
    try:
        op4, new = c.load('egg', alloc=5)
        assert False
    except TypeError:
        pass
    # after invalidation
    op1.set_element(1, 2, 5.2)
    c.invalidate_all()
    assert op1._array[1,2] == 0.0
    try:
        op4 = c.load('egg')
        assert False
    except KeyError:
        pass
    try:
        op4, new = c.load('egg', alloc=(lf, 'one_body', 5))
        assert False
    except TypeError:
        pass
    try:
        op4, new = c.load('egg', alloc=5)
        assert False
    except TypeError:
        pass
    op4, new = c.load('egg', alloc=(lf, 'one_body', 10))
    assert new
    assert op1 is op4
    op5 = c.load('egg')
    assert op1 is op5


def test_cache_basic_exceptions():
    c = Cache()
    try:
        c.load('boo')
        assert False
    except KeyError:
        pass

    c.dump('bar', np.zeros(4, float))
    try:
        c.load('bar', alloc=5)
        assert False
    except TypeError:
        pass

    try:
        c.load()
        assert False
    except TypeError:
        pass

    try:
        c.load('foo', sadfj=4)
        assert False
    except TypeError:
        pass

    try:
        c.load('foo', alloc=3, sdasffd=0)
        assert False
    except TypeError:
        pass

    try:
        c.load('foo', alloc=3, default=0)
        assert False
    except TypeError:
        pass

    try:
        c.load('foo', jgfjg=3, default=0)
        assert False
    except TypeError:
        pass

    try:
        c.dump()
        assert False
    except TypeError:
        pass

    try:
        c.dump('one')
        assert False
    except TypeError:
        pass


def test_discard():
    c = Cache()
    c.dump('foo', 5)
    c.dump('bar', 6)
    c.discard('foo')
    assert not c.has('foo')
    assert c.has('bar')


def check_array_store_common(mode, fn):
    with ArrayStore.from_mode(mode, fn) as store:
        arr1 = np.random.normal(0, 1, (10, 10))
        arr2 = np.zeros((10, 10))
        store.dump(arr1, 'aaaa', 5)
        if store.load(arr2, 'aaaa', 5):
            assert (arr1 == arr2).all()
        arr3 = np.random.normal(0, 1, (10, 10))
        store.dump(arr3, 'aaaa', 5,)
        if store.load(arr2, 'aaaa', 5):
            assert (arr3 == arr2).all()
        store.dump(arr3, *('aaaa', 5))
        if store.load(arr2, 'aaaa', 5):
            assert (arr3 == arr2).all()
            assert ('aaaa', 5) in store
            store.rename(('aaaa', 5), ('b', '1'))
            assert ('b', '1') in store
            assert ('aaaa', 5) not in store
        arr4 = np.random.normal(0, 1, (10, 10))
        if store.load(arr4, 'b', 1):
            assert (arr3 == arr4).all()
            try:
                arr5 = np.random.normal(0, 1, (11, 11))
                store.dump(arr5, 'b', 1)
                assert False
            except TypeError:
                pass
        try:
            store.load('b', 1, foobar=None)
            assert False
        except TypeError:
            pass


def test_array_store_disk():
    dn = tempfile.mkdtemp('horton.test.test_cache.test_array_store_disk')
    fn = '%s/test.h5' % dn
    try:
        check_array_store_common('disk', fn)
        assert not os.path.isfile(fn)
    finally:
        if os.path.isfile(fn):
           os.remove(fn)
        os.rmdir(dn)


def test_array_store_core():
    fn = 'horton.test.test_cache.test_array_store_core'
    check_array_store_common('core', fn)
    assert not os.path.isfile(fn)

def test_array_store_fake():
    check_array_store_common('fake', None)
