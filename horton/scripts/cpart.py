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
'''Utility functions for the ``horton-cpart.py`` script'''


__all__ = ['cpart_schemes']


def get_cpart_schemes():
    '''Return a dictionary with all cpart schemes'''
    import horton.part
    from horton.part.base import CPart
    cpart_schemes = {}
    for o in vars(horton.part).itervalues():
        if isinstance(o, type) and issubclass(o, CPart) and o.name is not None:
            cpart_schemes[o.name] = o
    return cpart_schemes


cpart_schemes = get_cpart_schemes()
