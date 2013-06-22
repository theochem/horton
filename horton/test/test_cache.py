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


def test_basics1():
    c = Cache()
    c.dump('foo', 5)
    assert c['foo'] == 5
    assert c[('foo',)] == 5
    assert c.load('foo') == 5
    assert c.load(('foo',)) == 5
    c.dump('foo', 4, 6)
    assert c['foo', 4] == 6
    assert c[('foo', 4)] == 6
    assert c.load('foo', 4) == 6
    assert c.load(('foo', 4)) == 6
    assert len(c) == 2
    c.clear()
    assert len(c._store) == 0
    assert len(c) == 0


def test_basics2():
    c = Cache()
    c['foo'] = 5
    assert c['foo'] == 5
    assert c[('foo',)] == 5
    assert c.load('foo') == 5
    assert c.load(('foo',)) == 5
    c['foo', 4] = 6
    assert c['foo', 4] == 6
    assert c[('foo', 4)] == 6
    assert c.load('foo', 4) == 6
    assert c.load(('foo', 4)) == 6
    assert len(c) == 2
    c.clear()
    assert len(c._store) == 0
    assert len(c) == 0


def test_alloc1():
    c = Cache()
    assert 'bar' not in c
    tmp, new = c.load('bar', alloc=5)
    assert new
    assert (tmp == 0).all()
    assert tmp.shape == (5,)
    assert issubclass(tmp.dtype.type, float)
    assert 'bar' in c
    tmp[3] = 1
    bis = c.load('bar')
    assert bis is tmp
    assert bis[3] == 1
    tris, new = c.load('bar', alloc=5)
    assert not new
    assert tris is tmp


def test_alloc2():
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


def test_multiple():
    c = Cache()
    assert ('a', 1) not in c
    c.dump('a', 1, 2)
    assert ('a', 1) in c
    assert c.load('a', 1) == 2


def test_allocation():
    c = Cache()
    assert 'egg' not in c
    tmp, new = c.load('egg', alloc=(5,10))
    assert new
    assert (tmp == 0).all()
    assert tmp.shape == (5,10)
    assert issubclass(tmp.dtype.type, float)
    assert 'egg' in c
    assert 'bar' not in c
    tmp[:] = 1.0
    c.clear()
    assert 'egg' not in c
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
    assert 'egg' in c
    # simple load should now work
    tris = c.load('egg')
    assert tris is tmp


def test_default():
    c = Cache()
    # with scalars
    assert c.load('egg', default=5)
    c.dump('egg', 5)
    assert c.load('egg') == 5
    assert c.load('egg', default=6) == 5
    c.clear()
    assert c.load('egg', default=6) == 6
    try:
        c.load('egg')
        assert False
    except KeyError:
        pass
    c.clear()
    assert c.load('egg', default=None) == None
    try:
        c.load('egg')
        assert False
    except KeyError:
        pass
    # with arrays
    c.dump('floep', np.array([3.1, 5.1]))
    assert (c.load('floep', default=3) == np.array([3.1, 5.1])).all()
    c.clear()
    assert c.load('floep', default=3) == 3
    try:
        c.load('floep')
        assert False
    except KeyError:
        pass


def test_dense_expansion():
    from horton.matrix import DenseLinalgFactory, DenseExpansion
    lf = DenseLinalgFactory()
    c = Cache()
    op1, new = c.load('egg', alloc=(lf, 'expansion', 10, 9))
    assert new
    assert isinstance(op1, DenseExpansion)
    assert op1.nbasis == 10
    op2 = c.load('egg')
    assert op1 is op2
    op3, new = c.load('egg', alloc=(lf, 'expansion', 10, 9))
    assert not new
    assert op1 is op3
    # things that should not work
    try:
        op4, new = c.load('egg', alloc=(lf, 'expansion', 5))
        assert False
    except TypeError:
        pass
    try:
        op4, new = c.load('egg', alloc=(lf, 'expansion', 10, 5))
        assert False
    except TypeError:
        pass
    try:
        op4, new = c.load('egg', alloc=5)
        assert False
    except TypeError:
        pass
    # after invalidation
    op1.coeffs[1, 2] = 5.2
    c.clear()
    assert op1.coeffs[1,2] == 0.0
    try:
        op4 = c.load('egg')
        assert False
    except KeyError:
        pass
    try:
        op4, new = c.load('egg', alloc=(lf, 'expansion', 5))
        assert False
    except TypeError:
        pass
    try:
        op4, new = c.load('egg', alloc=(lf, 'expansion', 10, 5))
        assert False
    except TypeError:
        pass
    try:
        op4, new = c.load('egg', alloc=5)
        assert False
    except TypeError:
        pass
    op4, new = c.load('egg', alloc=(lf, 'expansion', 10))
    assert new
    assert op1 is op4
    op5 = c.load('egg')
    assert op1 is op5


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
    c.clear()
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


def test_dense_two_body():
    from horton.matrix import DenseLinalgFactory, DenseTwoBody
    lf = DenseLinalgFactory()
    c = Cache()
    op1, new = c.load('egg', alloc=(lf, 'two_body', 10))
    assert new
    assert isinstance(op1, DenseTwoBody)
    assert op1.nbasis == 10
    op2 = c.load('egg')
    assert op1 is op2
    op3, new = c.load('egg', alloc=(lf, 'two_body', 10))
    assert not new
    assert op1 is op3
    # things that should not work
    try:
        op4, new = c.load('egg', alloc=(lf, 'two_body', 5))
        assert False
    except TypeError:
        pass
    try:
        op4, new = c.load('egg', alloc=5)
        assert False
    except TypeError:
        pass
    # after invalidation
    op1.set_element(1, 2, 1, 2, 5.2)
    c.clear()
    assert op1._array[1,2,1,2] == 0.0
    try:
        op4 = c.load('egg')
        assert False
    except KeyError:
        pass
    try:
        op4, new = c.load('egg', alloc=(lf, 'two_body', 5))
        assert False
    except TypeError:
        pass
    try:
        op4, new = c.load('egg', alloc=5)
        assert False
    except TypeError:
        pass
    op4, new = c.load('egg', alloc=(lf, 'two_body', 10))
    assert new
    assert op1 is op4
    op5 = c.load('egg')
    assert op1 is op5


def test_basic_exceptions():
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

    try:
        c.clear_item()
        assert False
    except TypeError:
        pass


def test_dealloc():
    c = Cache()
    c.dump('foo', 5)
    c.dump('bar', 6)
    c.clear_item('foo', dealloc=True)
    assert 'foo' not in c
    assert 'bar' in c
    assert len(c._store) == 1
    c.dump('foo', 5)
    c.clear(dealloc=True)
    assert 'foo' not in c
    assert 'bar' not in c
    assert len(c._store) == 0


def test_dump_unpack():
    c = Cache()
    c.dump(('foo',), 5)
    assert 'foo' in c


def test_iter():
    c = Cache()
    c.dump('foo', 5)
    c.dump('bar', 6)
    assert sorted(c.iterkeys()) == ['bar', 'foo']
    assert sorted(c.itervalues()) == [5, 6]
    assert sorted(c.iteritems()) == [('bar', 6), ('foo', 5)]
    assert len(c) == 2
    assert sorted(c) == ['bar', 'foo']
