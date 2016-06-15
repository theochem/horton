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
'''Physicochemical constants in atomic units

   These are the physical constants defined in this module (in atomic units):

'''

boltzmann = 3.1668154051341965e-06
avogadro = 6.0221415e23
lightspeed = 137.03599975303575
planck = 6.2831853071795864769


# automatically spice up the docstrings

lines = [
    '    ================  ==================',
    '    Name              Value             ',
    '    ================  ==================',
]

for key, value in sorted(globals().iteritems()):
    if not isinstance(value, float):
        continue
    lines.append('    %16s  %.10e' % (key, value))
lines.append('    ================  ==================')

__doc__ += '\n'.join(lines)
