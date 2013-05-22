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


from horton.test.common import check_script


def test_example_001_hf_water():
    check_script('./run.py', 'examples/001_hf_water')


def test_example_002_hfs_water():
    check_script('./run.py', 'examples/002_hfs_water')


def test_example_003_o3lyp_water():
    check_script('./run.py', 'examples/003_o3lyp_water')


def test_example_004_wpart():
    check_script('./run.py', 'examples/004_wpart')
