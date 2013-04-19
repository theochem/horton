# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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
'''Package for density-based partitioning (fuzzy atoms)'''

from horton.part.base import *
from horton.part.becke import *
from horton.part.hirshfeld import *
from horton.part.hirshfeld_i import *
from horton.part.hirshfeld_e import *
from horton.part.linalg import *
from horton.part.proatomdb import *
from horton.part.stockholder import *

wpart_schemes = {}
for o in globals().values():
    if isinstance(o, type) and issubclass(o, WPart) and o.name is not None:
        wpart_schemes[o.name] = o

cpart_schemes = {}
for o in globals().values():
    if isinstance(o, type) and issubclass(o, CPart) and o.name is not None:
        cpart_schemes[o.name] = o

del o
