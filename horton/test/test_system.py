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


def test_load_water_number():
    fn = context.get_fn('test/water_number.xyz')
    system = System.from_file(fn)
    check_water(system)


def test_load_water_element():
    fn = context.get_fn('test/water_element.xyz')
    system = System.from_file(fn)
    check_water(system)


def check_water(system):
    assert system.numbers[0] == 1
    assert system.numbers[1] == 8
    assert system.numbers[2] == 1
    assert abs(np.linalg.norm(system.coordinates[0] - system.coordinates[1])/angstrom - 0.96) < 1e-5
    assert abs(np.linalg.norm(system.coordinates[2] - system.coordinates[1])/angstrom - 0.96) < 1e-5
