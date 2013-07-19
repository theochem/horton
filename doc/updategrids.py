#!/usr/bin/env python
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


from glob import glob
from cStringIO import StringIO
import os
from horton import periodic

from common import write_if_changed

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

f = StringIO()

def print_markers():
    print >> f, '======= === ====', ' '.join(['=========']*len(grids))

print_markers()
print >> f, 'Element   Z Zeff', ' '.join([grid.ljust(9) for grid in grids])
print >> f, '\         \ \   ', ' '.join([alt.ljust(9) for alt in alts])
print_markers()
for z, zeff in elements:
    print >> f, '%7s %3i %4i' % (periodic[z].symbol, z, zeff),
    for grid in grids:
        npoint = npoints.get((z, zeff, grid))
        if npoint is None:
            print >> f, '--       ',
        else:
            print >> f, '%9i' % npoint,
    print >> f
print_markers()

write_if_changed('grids.rst.inc', f.getvalue())
