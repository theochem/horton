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


import numpy as np
from nose.tools import assert_raises

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import


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
    e.clear()
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
    ar1, new = c.load('egg', alloc=(5,10))
    assert new
    assert (ar1 == 0).all()
    assert ar1.shape == (5,10)
    assert issubclass(ar1.dtype.type, float)
    assert 'egg' in c
    assert 'bar' not in c
    with assert_raises(TypeError):
        c.load('egg', alloc=10)
    with assert_raises(TypeError):
        c.load('egg', alloc=(10,5))
    ar1[:] = 1.0
    c.clear()
    assert 'egg' not in c
    assert (ar1[:] == 0.0).all()
    # try to load it, while it is no longer valid
    with assert_raises(KeyError):
        ar2 = c.load('egg')
    # properly load it anew
    ar2, new = c.load('egg', alloc=(5,10))
    assert new
    assert ar2 is ar1 # still the same array, just cleared.
    assert 'egg' in c
    # simple load should now work
    ar3 = c.load('egg')
    assert ar3 is ar1
    # clear again and use different alloc
    c.clear()
    ar4, new = c.load('egg', alloc=(5,1,2))
    assert new
    assert ar4.shape == (5,1,2)
    assert not ar4 is ar1


def test_default():
    c = Cache()
    # with scalars
    assert c.load('egg', default=5)
    c.dump('egg', 5)
    assert c.load('egg') == 5
    assert c.load('egg', default=6) == 5
    c.clear()
    assert c.load('egg', default=6) == 6
    with assert_raises(KeyError):
        c.load('egg')
    c.clear()
    assert c.load('egg', default=None) == None
    with assert_raises(KeyError):
        c.load('egg')
    # with arrays
    c.dump('floep', np.array([3.1, 5.1]))
    assert (c.load('floep', default=3) == np.array([3.1, 5.1])).all()
    c.clear()
    assert c.load('floep', default=3) == 3
    with assert_raises(KeyError):
        c.load('floep')


def test_dense_expansion():
    from horton.matrix import DenseLinalgFactory, DenseExpansion
    lf = DenseLinalgFactory()
    c = Cache()
    exp1, new = c.load('egg', alloc=(lf.create_expansion, 10, 9))
    assert new
    assert isinstance(exp1, DenseExpansion)
    assert exp1.nbasis == 10
    assert exp1.nfn == 9
    exp2 = c.load('egg')
    assert exp1 is exp2
    exp3, new = c.load('egg', alloc=(lf.create_expansion, 10, 9))
    assert not new
    assert exp1 is exp3
    # things that should not work
    with assert_raises(TypeError):
        exp4, new = c.load('egg', alloc=(lf.create_expansion, 5))
    with assert_raises(TypeError):
        exp4, new = c.load('egg', alloc=(lf.create_expansion, 10, 5))
    with assert_raises(TypeError):
        exp4, new = c.load('egg', alloc=5)
    # after clearing
    exp1.coeffs[1, 2] = 5.2
    c.clear()
    assert exp1.coeffs[1,2] == 0.0
    with assert_raises(KeyError):
        exp4 = c.load('egg')
    exp4, new = c.load('egg', alloc=(lf.create_expansion, 10, 9))
    assert new
    assert exp1 is exp4
    exp5 = c.load('egg')
    assert exp1 is exp5
    # default_nbasis
    lf.default_nbasis = 5
    with assert_raises(TypeError):
        c.load('egg', alloc=lf.create_expansion)
    c.clear()
    exp6, new = c.load('egg', alloc=lf.create_expansion)
    assert new
    assert not exp5 is exp6
    assert exp6.nbasis == 5
    assert exp6.nfn == 5


def test_dense_two_index():
    from horton.matrix import DenseLinalgFactory, DenseTwoIndex
    lf = DenseLinalgFactory()
    c = Cache()
    op1, new = c.load('egg', alloc=(lf.create_two_index, 10))
    assert new
    assert isinstance(op1, DenseTwoIndex)
    assert op1.nbasis == 10
    op2 = c.load('egg')
    assert op1 is op2
    op3, new = c.load('egg', alloc=(lf.create_two_index, 10))
    assert not new
    assert op1 is op3
    # things that should not work
    with assert_raises(TypeError):
        op4, new = c.load('egg', alloc=(lf.create_two_index, 5))
    with assert_raises(TypeError):
        op4, new = c.load('egg', alloc=5)
    # after clearing
    op1.set_element(1, 2, 5.2)
    c.clear()
    assert op1._array[1,2] == 0.0
    with assert_raises(KeyError):
        op4 = c.load('egg')
    op4, new = c.load('egg', alloc=(lf.create_two_index, 10))
    assert new
    assert op1 is op4
    op5 = c.load('egg')
    assert op1 is op5
    # default_nbasis
    lf.default_nbasis = 5
    with assert_raises(TypeError):
        c.load('egg', alloc=lf.create_two_index)
    c.clear()
    op6, new = c.load('egg', alloc=lf.create_two_index)
    assert new
    assert not op5 is op6
    assert op6.nbasis == 5
    # the new method of the two-index object
    op7, new = c.load('bork', alloc=op6.new)
    assert new
    assert not op5 is op7
    assert op7.nbasis == 5


def test_dense_four_index():
    from horton.matrix import DenseLinalgFactory, DenseFourIndex
    lf = DenseLinalgFactory()
    c = Cache()
    op1, new = c.load('egg', alloc=(lf.create_four_index, 10))
    assert new
    assert isinstance(op1, DenseFourIndex)
    assert op1.nbasis == 10
    op2 = c.load('egg')
    assert op1 is op2
    op3, new = c.load('egg', alloc=(lf.create_four_index, 10))
    assert not new
    assert op1 is op3
    # things that should not work
    with assert_raises(TypeError):
        op4, new = c.load('egg', alloc=(lf.create_four_index, 5))
    with assert_raises(TypeError):
        op4, new = c.load('egg', alloc=5)
    # after clearing
    op1.set_element(1, 2, 1, 2, 5.2)
    c.clear()
    assert op1._array[1,2,1,2] == 0.0
    with assert_raises(KeyError):
        op4 = c.load('egg')
    op4, new = c.load('egg', alloc=(lf.create_four_index, 10))
    assert new
    assert op1 is op4
    op5 = c.load('egg')
    assert op1 is op5
    # default_nbasis
    lf.default_nbasis = 5
    with assert_raises(TypeError):
        c.load('egg', alloc=lf.create_four_index)
    c.clear()
    op6, new = c.load('egg', alloc=lf.create_two_index)
    assert new
    assert not op5 is op6
    assert op6.nbasis == 5


def test_basic_exceptions():
    c = Cache()
    with assert_raises(KeyError):
        c.load('boo')
    c.dump('bar', np.zeros(4, float))
    with assert_raises(TypeError):
        c.load('bar', alloc=5)
    with assert_raises(TypeError):
        c.load()
    with assert_raises(TypeError):
        c.load('foo', sadfj=4)
    with assert_raises(TypeError):
        c.load('foo', alloc=3, sdasffd=0)
    with assert_raises(TypeError):
        c.load('foo', alloc=3, default=0)
    with assert_raises(TypeError):
        c.load('foo', jgfjg=3, default=0)
    with assert_raises(TypeError):
        c.dump()
    with assert_raises(TypeError):
        c.dump('one')
    with assert_raises(TypeError):
        c.clear_item()


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


def test_iter_tags():
    c = Cache()
    c.dump('foo', 5, tags='c')
    c.dump('bar', 6)
    c.dump('egg', 7, tags='op')
    c.dump('spam', 8, tags='co')
    assert sorted(c.iterkeys()) == ['bar', 'egg', 'foo', 'spam']
    assert sorted(c.itervalues()) == [5, 6, 7, 8]
    assert sorted(c.iteritems()) == [('bar', 6), ('egg', 7), ('foo', 5), ('spam', 8)]
    assert sorted(c.iterkeys(tags='c')) == ['foo', 'spam']
    assert sorted(c.itervalues(tags='c')) == [5, 8]
    assert sorted(c.iteritems(tags='c')) == [('foo', 5), ('spam', 8)]
    assert sorted(c.iterkeys(tags='a')) == []
    assert sorted(c.itervalues(tags='a')) == []
    assert sorted(c.iteritems(tags='a')) == []
    assert len(c) == 4
    assert sorted(c) == ['bar', 'egg', 'foo', 'spam']


def test_tags():
    c = Cache()
    c.dump('a', 5, tags='ab')
    # In a normal load call, the tags should not be allowed
    with assert_raises(TypeError):
        assert c.load('a', tags='a') == 5
    with assert_raises(TypeError):
        assert c.load('a', tags='ab') == 5
    with assert_raises(TypeError):
        assert c.load('a', tags='abc') == 5
    with assert_raises(TypeError):
        assert c.load('b', default=5, tags='abc') == 5
    # clear with other tags
    c.clear(tags='cd')
    assert len(c) == 1
    # clear with correct tag
    c.clear(tags='a')
    assert len(c) == 0
    # use in combination with alloc
    tmp1, new = c.load('tmp', alloc=5, tags='qw')
    assert new
    tmp2, new = c.load('tmp', alloc=5, tags='qw')
    assert not new
    assert tmp1 is tmp2
    with assert_raises(ValueError):
        c.load('tmp', alloc=5, tags='w')
    with assert_raises(ValueError):
        c.load('tmp', alloc=5, tags='aw')
    with assert_raises(ValueError):
        c.load('tmp', alloc=5, tags='ab')
    c.clear()
    tmp3, new = c.load('tmp', alloc=5, tags='qw')
    assert new
    assert tmp1 is tmp3
    with assert_raises(ValueError):
        c.load('tmp', alloc=5, tags='w')
    with assert_raises(ValueError):
        c.load('tmp', alloc=5, tags='aw')
    with assert_raises(ValueError):
        c.load('tmp', alloc=5, tags='ab')
