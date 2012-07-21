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
    olp1 = sys1.get_overlap()
    olp2 = sys2.get_overlap()
    mask = abs(olp1._array) > 1e-5
    delta = olp1._array - olp2._array
    expect = olp1._array
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
    kin1 = sys1.get_kinetic()
    sys2 = System.from_file(fn_fchk)
    kin2 = sys2.get_kinetic()
    mask = abs(kin1._array) > 1e-5
    delta = kin1._array - kin2._array
    expect = kin1._array
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
    na1 = sys1.get_nuclear_attraction()
    na2 = sys2.get_nuclear_attraction()
    mask = abs(na1._array) > 1e-5
    expect = na1._array
    result = na2._array
    delta = expect - result
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



def check_electron_repulsion(fn_fchk, check_zeros=False):
    fn_log = fn_fchk[:-5] + '.log'
    sys1 = System.from_file(fn_fchk, fn_log)
    er1 = sys1.get_electron_repulsion()
    sys2 = System.from_file(fn_fchk)
    er2 = sys2.get_electron_repulsion()
    mask = abs(er1._array) > 1e-6
    expect = er1._array
    got = er2._array
    if check_zeros:
        assert ((expect == 0.0) == (got == 0.0)).all()
    delta = expect - got
    error = (delta[mask]/expect[mask]).max()
    print error
    assert error < 1e-5


def test_electron_repulsion_water_sto3g_hf():
    check_electron_repulsion(context.get_fn('test/water_sto3g_hf_g03.fchk'), True)


def test_electron_repulsion_water_ccpvdz_pure_hf():
    check_electron_repulsion(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'))


def test_electron_repulsion_water_ccpvdz_cart_hf():
    check_electron_repulsion(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'))
