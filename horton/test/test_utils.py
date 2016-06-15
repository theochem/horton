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

from horton import *  # pylint: disable=wildcard-import,unused-wildcard-import

def test_typecheck_geo():
    numbers = np.array([11, 12, 13, 14])
    coordinates = np.random.normal(0, 1, (4, 3))
    pseudo_numbers = np.array([1, 2, 3, 4])

    result = typecheck_geo(coordinates, numbers)
    assert result[0] == 4
    assert (result[1] == coordinates).all()
    assert (result[2] == numbers).all()
    assert (result[3] == numbers).all()
    assert issubclass(result[3].dtype.type, float)

    result = typecheck_geo(pseudo_numbers=pseudo_numbers, need_coordinates=False, need_numbers=False)
    assert result != 4 # this was a bug
    assert result[0] == 4
    assert (result[1] == pseudo_numbers).all()
    assert issubclass(result[1].dtype.type, float)
