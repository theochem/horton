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


import os

from horton import System, context
from horton.scripts.wpart import *
from horton.part.test.common import get_proatomdb_ref


def test_parse_grid_1():
    for sgrid in 'coarse', 'medium', 'fine', 'veryfine':
        atspecs = parse_grid(sgrid, None, None)
        fn = context.get_fn('grids/%s.txt' % atspecs)
        assert os.path.isfile(fn)


def test_parse_grid_2():
    sys = System.from_file(context.get_fn('test/water_number.xyz'))
    padb = get_proatomdb_ref([1, 8], 1, 1)
    for sgrid in '6', '38', '110':
        nll = int(sgrid)
        atspecs = parse_grid(sgrid, sys, padb)
        for i in xrange(sys.natom):
            assert len(atspecs[i]) == 3
            assert atspecs[i][0] == padb.get_rtransform(sys.numbers[i])
            assert atspecs[i][2] == nll
