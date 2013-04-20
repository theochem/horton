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


from horton.scripts.atomdb import *


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
    try:
        list(iter_elements('8-6,2')) == [2]
        assert False
    except ValueError:
        pass


def test_iter_mults():
    assert list(iter_mults(4, True)) == [1]
    assert list(iter_mults(4, False)) == [1, 3]
    assert list(iter_mults(42, True)) == [7]
    assert list(iter_mults(42, False)) == [7, 5, 3, 1]
    assert list(iter_mults(47, True)) == [2]
    assert list(iter_mults(47, False)) == [2]


def test_iter_states():
    l = list(iter_states('1-3', 2, 2, True))
    assert l == [(1, -2, 2), (1, -1, 1), (1, 0, 2), (2, -2, 1), (2, -1, 2),
                 (2, 0, 1), (2, 1, 2), (3, -2, 2), (3, -1, 1), (3, 0, 2),
                 (3, 1, 1), (3, 2, 2)]
    l = list(iter_states('1, 6-8', 3, 1, False))
    assert l == [(1, -1, 1), (1, -1, 3), (1, 0, 2), (6, -1, 4), (6, -1, 2),
                 (6, 0, 3), (6, 0, 5), (6, 0, 1), (6, 1, 2), (6, 1, 4),
                 (6, 2, 1), (6, 2, 3), (6, 3, 2), (6, 3, 4), (7, -1, 3),
                 (7, -1, 1), (7, 0, 4), (7, 0, 2), (7, 1, 3), (7, 1, 5),
                 (7, 1, 1), (7, 2, 2), (7, 2, 4), (7, 3, 1), (7, 3, 3),
                 (8, -1, 2), (8, 0, 3), (8, 0, 1), (8, 1, 4), (8, 1, 2),
                 (8, 2, 3), (8, 2, 5), (8, 2, 1), (8, 3, 2), (8, 3, 4)]
