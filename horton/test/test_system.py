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


def test_overlap_water_sto3g_hf():
    sys1 = System.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'), context.get_fn('test/water_sto3g_hf_g03.log'))
    sys2 = System.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    sys2.init_overlap()
    error = abs(sys1.operators['olp']._array - sys2.operators['olp']._array).max()
    assert error < 1e-6


def test_overlap_water_ccpvdz_pure_hf():
    sys1 = System.from_file(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'), context.get_fn('test/water_ccpvdz_pure_hf_g03.log'))
    sys2 = System.from_file(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'))
    sys2.init_overlap()
    error = abs(sys1.operators['olp']._array - sys2.operators['olp']._array).max()
    assert error < 1e-6


def test_overlap_water_ccpvdz_cart_hf():
    sys1 = System.from_file(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'), context.get_fn('test/water_ccpvdz_cart_hf_g03.log'))
    sys2 = System.from_file(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'))
    sys2.init_overlap()
    error = abs(sys1.operators['olp']._array - sys2.operators['olp']._array).max()
    assert error < 1e-6


def test_overlap_co_ccpv5z_pure_hf():
    sys1 = System.from_file(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'), context.get_fn('test/co_ccpv5z_pure_hf_g03.log'))
    sys2 = System.from_file(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'))
    sys2.init_overlap()
    delta = sys1.operators['olp']._array - sys2.operators['olp']._array
    for index, row in enumerate((abs(delta) > 1e-4).astype(int)):
        print index, ''.join({0:' ',1:'#'}[i] for i in row)
    print 'Gaussian'
    print sys1.operators['olp']._array[80,:10]
    print 'Horton'
    print sys2.operators['olp']._array[80,:10]
    error = abs(sys1.operators['olp']._array - sys2.operators['olp']._array).max()
    assert error < 1e-6


def test_overlap_co_ccpv5z_cart_hf():
    sys1 = System.from_file(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'), context.get_fn('test/co_ccpv5z_cart_hf_g03.log'))
    sys2 = System.from_file(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'))
    sys2.init_overlap()
    delta = sys1.operators['olp']._array - sys2.operators['olp']._array
    error = abs(sys1.operators['olp']._array - sys2.operators['olp']._array).max()
    assert error < 1e-6
