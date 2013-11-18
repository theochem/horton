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


import numpy as np

from horton import *


def test_cell():
    fn = context.get_fn('test/water_element.xyz')
    system = System.from_file(fn, cell=Cell(np.identity(3, float)*10))
    assert abs(system.cell.rvecs - np.identity(3, float)*10).max() < 1e-10


def test_nucnuc():
    fn_fchk = context.get_fn('test/hf_sto3g.fchk')
    sys = System.from_file(fn_fchk)
    assert abs(sys.compute_nucnuc() - 4.7247965053) < 1e-5


def check_scf_dms(fn_fchk):
    sys = System.from_file(fn_fchk)
    assert abs(sys.wfn.dm_full._array - sys.cache['dm_scf_full']._array).max() < 1e-7
    if 'dm_scf_spin' in sys.cache:
        assert abs(sys.wfn.dm_spin._array - sys.cache['dm_scf_spin']._array).max() < 1e-7


def test_scf_dms_water_sto3g_hf():
    check_scf_dms(context.get_fn('test/water_sto3g_hf_g03.fchk'))

def test_scf_dms_lih_321g_hf_TT():
    check_scf_dms(context.get_fn('test/li_h_3-21G_hf_g09.fchk'))

def test_scf_dms_water_ccpvdz_pure_hf():
    check_scf_dms(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'))

def test_scf_dms_water_ccpvdz_cart_hf():
    check_scf_dms(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'))


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
    na1 = sys1.get_nuclear_attraction()
    na2 = sys2.get_nuclear_attraction()
    mask = abs(na1._array) > 1e-5
    expect = na1._array
    result = na2._array
    delta = -expect - result
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
    assert error < 1e-5


def test_electron_repulsion_water_sto3g_hf():
    check_electron_repulsion(context.get_fn('test/water_sto3g_hf_g03.fchk'), True)


def test_electron_repulsion_water_ccpvdz_pure_hf():
    check_electron_repulsion(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'))


def test_electron_repulsion_water_ccpvdz_cart_hf():
    check_electron_repulsion(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'))


def check_grid_fn(fn_fchk, use_output_arg):
    sys = System.from_file(fn_fchk)
    grid = BeckeMolGrid(sys, 'tv-13.7-4', random_rotate=False)
    if use_output_arg:
        rhos = np.zeros(grid.size)
        sys.compute_grid_density(grid.points, rhos)
    else:
        rhos = sys.compute_grid_density(grid.points)
    pop = grid.integrate(rhos)
    assert abs(pop-sys.wfn.nel) < 2e-3


def test_grid_fn_h_sto3g_T():
    check_grid_fn(context.get_fn('test/h_sto3g.fchk'), True)

def test_grid_fn_h_sto3g_F():
    check_grid_fn(context.get_fn('test/h_sto3g.fchk'), False)


def test_grid_fn_lih_321g_hf_T():
    check_grid_fn(context.get_fn('test/li_h_3-21G_hf_g09.fchk'), True)


def test_grid_fn_water_sto3g_hf_T():
    check_grid_fn(context.get_fn('test/water_sto3g_hf_g03.fchk'), True)


def test_grid_fn_water_ccpvdz_pure_hf_F():
    check_grid_fn(context.get_fn('test/water_ccpvdz_pure_hf_g03.fchk'), False)


def test_grid_fn_water_ccpvdz_cart_hf_F():
    check_grid_fn(context.get_fn('test/water_ccpvdz_cart_hf_g03.fchk'), False)


def test_grid_fn_co_ccpv5z_pure_hf_T():
    check_grid_fn(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'), True)

def test_grid_fn_co_ccpv5z_pure_hf_F():
    check_grid_fn(context.get_fn('test/co_ccpv5z_pure_hf_g03.fchk'), False)

def test_grid_fn_co_ccpv5z_cart_hf_T():
    check_grid_fn(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'), True)

def test_grid_fn_co_ccpv5z_cart_hf_F():
    check_grid_fn(context.get_fn('test/co_ccpv5z_cart_hf_g03.fchk'), False)


def check_normalization_dm_full_azirine(fn_fchk):
    system = System.from_file(fn_fchk)
    olp = system.get_overlap()
    dm = system.wfn.dm_full
    check_dm(dm, olp, system.lf, 'full', eps=1e-2, occ_max=2)
    assert abs(olp.expectation_value(dm) - 22.0) < 1e-3


def test_normalization_dm_full_azirine_ccd():
    check_normalization_dm_full_azirine(context.get_fn('test/2h-azirine-ccd.fchk'))


def test_normalization_dm_full_azirine_cis():
    check_normalization_dm_full_azirine(context.get_fn('test/2h-azirine-cis.fchk'))


def test_normalization_dm_full_azirine_mp2():
    check_normalization_dm_full_azirine(context.get_fn('test/2h-azirine-mp2.fchk'))


def test_normalization_dm_full_azirine_mp3():
    check_normalization_dm_full_azirine(context.get_fn('test/2h-azirine-mp3.fchk'))


def test_update_obasis():
    sys = System.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    assert sys.obasis_desc is None
    assert 'energy' in sys.extra
    olp = sys.get_overlap()
    assert olp.nbasis == 7
    sys.update_obasis('3-21G')
    assert sys.obasis_desc.default == '3-21G'
    assert 'olp' not in sys.cache
    assert len(sys.extra) == 0
    olp = sys.get_overlap()
    assert olp.nbasis == 13


def test_update_coordinates1():
    sys = System.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    olp = sys.get_overlap()
    assert len(sys.extra) > 0
    assert len(sys.cache) > 0
    sys.update_coordinates(np.array([[0.0, 0.1, 0.2], [0.3, 0.4, 0.5], [0.6, 0.7, 0.8]]))
    assert sys.coordinates[1,2] == 0.5
    assert len(sys.extra) == 0
    assert len(sys.cache) == 0
    assert sys.obasis.centers[1,2] == 0.5


def test_update_coordinates2():
    sys = System.from_file(context.get_fn('test/water_sto3g_hf_g03.fchk'))
    olp = sys.get_overlap()
    assert len(sys.extra) > 0
    assert len(sys.cache) > 0
    sys.coordinates[:] = [[0.0, 0.1, 0.2], [0.3, 0.4, 0.5], [0.6, 0.7, 0.8]]
    sys.update_coordinates()
    assert sys.coordinates[1,2] == 0.5
    assert len(sys.extra) == 0
    assert len(sys.cache) == 0
    assert sys.obasis.centers[1,2] == 0.5
