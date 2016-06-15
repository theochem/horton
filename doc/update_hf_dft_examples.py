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
'''Update the table with HF/DFT examples based on files in data/examples/hf_dft.'''

from glob import glob
from cStringIO import StringIO
import json, os
from common import write_rst_table, write_if_changed

cases = []
for fn_py in sorted(glob('../data/examples/hf_dft/*.py')):
    with open(fn_py) as f:
        s = ''
        for line in f:
            if line.startswith('#JSON'):
                s += line[5:]
    if len(s) == 0:
        raise RuntimeError('Could not find JSON line in HF/DFT example script.')
    meta = json.loads(s)
    case = [meta['difficulty'], os.path.basename(fn_py), meta['lot'],
            meta['scf'], meta['description']]
    cases.append(case)

cases.sort()
table = [['File name', 'LOT', 'SCF', 'Description']] + [case[1:] for case in cases]

f = StringIO()
write_rst_table(f, table)
write_if_changed('hf_dft_examples.rst.inc', f.getvalue())
