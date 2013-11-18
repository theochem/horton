#!/usr/bin/env python
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
#pylint: skip-file


from nose.tools import assert_raises
from horton import *


def test_base_exceptions():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, mode='discard')
    with assert_raises(ValueError):
        # the default setting is local=true, which is not compatible with mode='discard'.
        dp = WPart(sys, grid)

    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, mode='keep')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(sys, grid)

    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, mode='only')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(sys, grid)

    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, mode='discard')
    with assert_raises(NotImplementedError):
        # It should not be possible to create instances of the base class.
        dp = WPart(sys, grid, local=False)


def test_wpart_schemes():
    assert 'b' in wpart_schemes
    assert 'h' in wpart_schemes
    assert 'hi' in wpart_schemes
    assert 'he' in wpart_schemes
    assert wpart_schemes['hi'] is HirshfeldIWPart
    assert wpart_schemes['hi'].options == ['lmax', 'threshold', 'maxiter', 'greedy', 'epsilon']
    assert not wpart_schemes['hi'].linear
    assert wpart_schemes['h'].linear
    assert wpart_schemes['b'].linear

    for WPartClass in wpart_schemes.itervalues():
        assert hasattr(WPartClass, 'options')


def test_cpart_schemes():
    assert 'h' in cpart_schemes
    assert 'hi' in cpart_schemes
    assert 'he' in cpart_schemes
    assert cpart_schemes['hi'] is HirshfeldICPart

    for CPartClass in cpart_schemes.itervalues():
        print CPartClass
        assert hasattr(CPartClass, 'options')
