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


import numpy as np
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
    c.invalidate()
    assert len(c._store) == 0


def test_cache_alloc1():
    c = Cache()
    tmp, new = c.load('bar', alloc=5)
    assert new
    assert (tmp == 0).all()
    assert tmp.shape == (5,)
    assert issubclass(tmp.dtype.type, float)
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
    tmp, new = c.load('egg', alloc=(5,10))
    assert new
    assert (tmp == 0).all()
    assert tmp.shape == (5,10)
    assert issubclass(tmp.dtype.type, float)
    tmp[:] = 1.0
    c.invalidate()
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
    # simple load should now work
    tris = c.load('egg')
    assert tris is tmp


def test_cache_exceptions():
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
        c.dump()
        assert False
    except TypeError:
        pass

    try:
        c.dump('one')
        assert False
    except TypeError:
        pass
