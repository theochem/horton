#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --


from glob import glob
from cStringIO import StringIO
import os
from horton import periodic

from common import write_rst_table, write_if_changed

elements = set([])
grids = set([])
npoints = {}

for fn_txt in glob('../data/grids/*.txt'):
    grid = os.path.basename(fn_txt)[:-4]
    grids.add(grid)
    with open(fn_txt) as f:
        for line in f:
            if line.startswith('# Number of points:'):
                npoint = int(line.split()[-1])
            words = line.split()
            if line[0] != '#' and len(words) == 2:
                assert isinstance(npoint, int)
                npoints[int(words[0]), int(words[1]), grid] = npoint
                npoint = None
                elements.add((int(words[0]), int(words[1])))


grids = sorted(grids)
alts = ['coarse', 'medium', 'fine', 'veryfine', 'ultrafine', 'insane']
elements = sorted(elements)

table = [
    ['Element', 'Z', 'Zeff'] + grids,
    [' ', ' ', ' '] + alts,
]
for z, zeff in elements:
    row = [periodic[z].symbol, str(z), str(zeff)]
    for grid in grids:
        npoint = npoints.get((z, zeff, grid))
        if npoint is None:
            row.append('--')
        else:
            row.append(str(npoint))
    table.append(row)

f = StringIO()
write_rst_table(f, table, nhead=2)
write_if_changed('grids.rst.inc', f.getvalue())
