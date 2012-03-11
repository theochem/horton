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


from horton import *


def test_fact():
    assert fact(-20) == 1
    assert fact(0) == 1
    assert fact(1) == 1
    assert fact(2) == 2
    assert fact(3) == 6
    assert fact(4) == 24
    assert fact(5) == 120


def test_fact2():
    assert fact2(-20) == 1
    assert fact2(0) == 1
    assert fact2(1) == 1
    assert fact2(2) == 2
    assert fact2(3) == 3
    assert fact2(4) == 8
    assert fact2(5) == 15
