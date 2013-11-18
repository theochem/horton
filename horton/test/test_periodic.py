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
    assert periodic[2].bs_radius is None
    assert periodic[' Be '].bs_radius == 1.05*angstrom
    assert periodic[' He '].wc_radius == 0.291*angstrom
    assert periodic[' C '].vdw_radius == 1.75*angstrom
