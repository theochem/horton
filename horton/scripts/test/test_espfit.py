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


from horton import angstrom
from horton.scripts.espfit import *

def test_wdens():
    assert parse_wdens('fubar.cube') == ('fubar.cube', 2e-4, 1.0)
    assert parse_wdens('fubar.cube:2e-3') == ('fubar.cube', 2e-3, 1.0)
    assert parse_wdens('fubar.cube:2e-2:0.5') == ('fubar.cube', 2e-2, 0.5)

def test_wnear():
    assert parse_wnear('1:1.0') == {1: (1.0*angstrom, 0.5*angstrom)}
    assert parse_wnear('1:1.0:0.3') == {1: (1.0*angstrom, 0.3*angstrom)}
    assert parse_wnear(['1:1.0', '2:1.2']) == {1: (1.0*angstrom, 0.5*angstrom), 2: (1.2*angstrom, 0.6*angstrom)}
    assert parse_wnear(['1:1.0:0.3', '2:1.2:0.2']) == {1: (1.0*angstrom, 0.3*angstrom), 2: (1.2*angstrom, 0.2*angstrom)}

def test_wfar():
    assert parse_wfar('4.3') == (4.3*angstrom, 1.0*angstrom)
    assert parse_wfar('4.2:0.3') == (4.2*angstrom, 0.3*angstrom)
