# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2015 The HORTON Development Team
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
#--
'''Utility functions for the ``horton-cpart.py`` script'''


import horton.part
from horton.part.base import CPart


__all__ = ['cpart_schemes']


def get_cpart_schemes():
    '''Return a dictionary with all cpart schemes'''
    cpart_schemes = {}
    for o in vars(horton.part).itervalues():
        if isinstance(o, type) and issubclass(o, CPart) and o.name is not None:
            cpart_schemes[o.name] = o
    return cpart_schemes


cpart_schemes = get_cpart_schemes()
