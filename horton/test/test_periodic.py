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


from horton import periodic
from horton.units import angstrom

def test_periodic():
    assert periodic['si'].number == 14
    assert periodic['He'].number == 2
    assert periodic['h'].symbol == 'H'
    assert periodic[3].symbol == 'Li'
    assert periodic['5'].symbol == 'B'
    assert periodic[' 5'].symbol == 'B'
    assert periodic[' B '].symbol == 'B'
    assert periodic[97].cov_radius is None
    assert periodic['cu'].cov_radius == 1.32*angstrom
    assert periodic['cu'].cov_radius_cordero == 1.32*angstrom
    assert periodic[2].cov_radius_bragg is None
    assert periodic[' Be '].cov_radius_bragg == 1.15*angstrom
    assert periodic[2].cov_radius_slater is None
    assert periodic[' Be '].cov_radius_slater == 1.05*angstrom
    assert periodic[' C '].vdw_radius_bondi == 1.7*angstrom
    assert periodic['\tAL'].vdw_radius_bondi is None
    assert periodic['Al'].vdw_radius_truhlar == 1.84*angstrom
    assert periodic['H'].vdw_radius_truhlar is None
    assert periodic['Al'].vdw_radius_rt == 2.03*angstrom
    assert periodic['100'].vdw_radius_rt is None
    assert periodic['nb'].vdw_radius_batsanov == 2.15*angstrom
    assert periodic['KR'].vdw_radius_batsanov is None
    assert periodic['cl'].vdw_radius_dreiding == 3.9503*angstrom/2
    assert periodic[18].vdw_radius_dreiding is None
    assert periodic['mG'].vdw_radius_uff == 3.021*angstrom/2
    assert periodic[118].vdw_radius_uff is None
    assert periodic[90].vdw_radius_mm3 == 2.74*angstrom
    assert periodic[100].vdw_radius_mm3 is None
    assert periodic[' He '].wc_radius == 0.291*angstrom
    assert periodic['Db'].wc_radius is None
    assert periodic['K'].cr_radius == 2.43*angstrom
    assert periodic['Fr'].cr_radius is None
    # check derived vdw radius
    assert periodic[1].vdw_radius == periodic[1].vdw_radius_bondi
    assert periodic[13].vdw_radius == periodic[13].vdw_radius_truhlar
    assert periodic['Sc'].vdw_radius == periodic['Sc'].vdw_radius_batsanov
    assert periodic[58].vdw_radius == periodic[58].vdw_radius_mm3
    # check derived becke radius
    assert periodic[1].becke_radius == periodic[1].cov_radius_slater
    assert periodic[2].becke_radius == periodic[2].cov_radius_cordero
