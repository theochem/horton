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


def test_get_max_nbasis():
    num_contractions = np.array([1, 2, 3, 1])
    con_types = np.array([2, 1, 0, -2, 3, 0, 1])
    assert get_max_nbasis(num_contractions, con_types) == 10

    num_contractions = np.array([1, 2, 3, 1])
    con_types = np.array([1, 1, 0, -2, -2, 0, 1])
    assert get_max_nbasis(num_contractions, con_types) == 6

    num_contractions = np.array([0, 1, 1])
    con_types = np.array([1, 1])
    try:
        get_max_nbasis(num_contractions, con_types)
        assert False
    except ValueError:
        pass

    num_contractions = np.array([1, 1, 1])
    con_types = np.array([-1, 1, 2])
    try:
        get_max_nbasis(num_contractions, con_types)
        assert False
    except ValueError:
        pass

    num_contractions = np.array([1, 1, 1])
    con_types = np.array([1, 1])
    try:
        get_max_nbasis(num_contractions, con_types)
        assert False
    except ValueError:
        pass
