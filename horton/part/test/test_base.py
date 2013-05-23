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


from horton import *


def test_base_exceptions():
    fn_fchk = context.get_fn('test/n2_hfs_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    rtf = ExpRTransform(1e-3, 1e1, 100)
    rgrid = RadialGrid(rtf)
    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, keep_subgrids=False)

    try:
        dp = WPart(sys, grid)
        assert False
    except ValueError:
        pass

    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, keep_subgrids=True)

    try:
        dp = WPart(sys, grid)
        assert False
    except NotImplementedError:
        pass

    grid = BeckeMolGrid(sys, (rgrid, 110), random_rotate=False, keep_subgrids=False)
    try:
        dp = WPart(sys, grid, local=False)
        assert False
    except NotImplementedError:
        pass


def test_wpart_schemes():
    assert 'b' in wpart_schemes
    assert 'h' in wpart_schemes
    assert 'hi' in wpart_schemes
    assert 'he' in wpart_schemes
    assert wpart_schemes['hi'] is HirshfeldIWPart
    assert wpart_schemes['hi'].options == ['threshold', 'maxiter', 'greedy']

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
