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
from horton.periodic import periodic
from horton.log import log
from horton.gbasis.gobasis import go_basis_families_list

from common import write_rst_table, write_if_changed


log.set_level(log.silent)


table = [['Basis set name', 'Supported elements']]
for family in go_basis_families_list:
    family.load()
    numbers = sorted(family.basis_atom_map.keys())
    s = []
    for i in xrange(len(numbers)):
        if (i > 0 and i < len(numbers)-1) and numbers[i-1] == numbers[i]-1 \
           and numbers[i+1] == numbers[i]+1:
            if s[-1] != '-':
                s.append('-')
        else:
            s.append(periodic[numbers[i]].symbol)
    range_str = ' '.join(s).replace(' - ', '--').replace(' ', ', ')
    table.append([family.name, range_str])

f = StringIO()
write_rst_table(f, table, nhead=1)
write_if_changed('basis.rst.inc', f.getvalue())
