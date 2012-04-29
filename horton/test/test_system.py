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



def check_overlap(fn_fchk):
    fn_log = fn_fchk[:-5] + '.log'
    sys1 = System.from_file(fn_fchk, fn_log)
    sys2 = System.from_file(fn_fchk)
    sys2.init_overlap()
    mask = abs(sys1.operators['olp']._array) > 1e-5
    delta = sys1.operators['olp']._array - sys2.operators['olp']._array
    expect = sys1.operators['olp']._array
    error = (delta[mask]/expect[mask]).max()
    assert error < 1e-5


def test_overlap_water_sto3g_hf():
    check_overlap(context.get_fn('test/water_sto3g_hf_g03.fchk'))


def test_overlap_water_ccpvdz_pure_hf():
    check_overlap(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'))


def test_overlap_water_ccpvdz_cart_hf():
    check_overlap(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'))


def test_overlap_co_ccpv5z_pure_hf():
    check_overlap(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'))


def test_overlap_co_ccpv5z_cart_hf():
    check_overlap(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'))




def check_kinetic(fn_fchk):
    fn_log = fn_fchk[:-5] + '.log'
    sys1 = System.from_file(fn_fchk, fn_log)
    sys2 = System.from_file(fn_fchk)
    sys2.init_kinetic()
    mask = abs(sys1.operators['kin']._array) > 1e-5
    delta = sys1.operators['kin']._array - sys2.operators['kin']._array
    expect = sys1.operators['kin']._array
    error = (delta[mask]/expect[mask]).max()
    assert error < 1e-5


def test_kinetic_water_sto3g_hf():
    check_kinetic(context.get_fn('test/water_sto3g_hf_g03.fchk'))


def test_kinetic_water_ccpvdz_pure_hf():
    check_kinetic(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'))


def test_kinetic_water_ccpvdz_cart_hf():
    check_kinetic(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'))


def test_kinetic_co_ccpv5z_pure_hf():
    check_kinetic(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'))


def test_kinetic_co_ccpv5z_cart_hf():
    check_kinetic(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'))


def check_nuclear_attraction(fn_fchk):
    fn_log = fn_fchk[:-5] + '.log'
    sys1 = System.from_file(fn_fchk, fn_log)
    sys2 = System.from_file(fn_fchk)
    charges = sys2.numbers.astype(float)
    centers = sys2.coordinates
    sys2.init_nuclear_attraction(charges, centers)
    mask = abs(sys1.operators['na']._array) > 1e-5
    expect = sys1.operators['na']._array
    result = sys2.operators['na']._array
    delta = expect - result
    np.set_printoptions(suppress=True)
    error = (delta[mask]/expect[mask]).max()
    assert error < 4e-5


def test_nuclear_attraction_water_sto3g_hf():
    check_nuclear_attraction(context.get_fn('test/water_sto3g_hf_g03.fchk'))


def test_nuclear_attraction_water_ccpvdz_pure_hf():
    check_nuclear_attraction(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'))


def test_nuclear_attraction_water_ccpvdz_cart_hf():
    check_nuclear_attraction(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'))


def test_nuclear_attraction_co_ccpv5z_pure_hf():
    check_nuclear_attraction(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'))


def test_nuclear_attraction_co_ccpv5z_cart_hf():
    check_nuclear_attraction(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'))
