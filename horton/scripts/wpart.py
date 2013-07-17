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
'''Code used by ``horton-wpart.py``'''


from horton.grid.cext import lebedev_laikov_npoints

__all__ = ['parse_grid']


def parse_grid(sgrid, system, proatomdb):
    if sgrid == 'coarse':
        return 'tv-13.7-3'
    elif sgrid == 'medium':
        return 'tv-13.7-4'
    elif sgrid == 'fine':
        return 'tv-13.7-5'
    elif sgrid == 'veryfine':
        return 'tv-13.7-6'
    elif sgrid == 'ultrafine':
        return 'tv-13.7-7'
    elif sgrid == 'insane':
        return 'tv-13.7-8'
    else:
        try:
            nll = int(sgrid)
        except ValueError:
            raise ValueError('If the grid argument is not coarse, medium, fine, veryfine, ultrafine or insane, it should be an integer.')

        if nll not in lebedev_laikov_npoints:
            raise ValueError('If the grid argument is an integer, it should be one of the following: %s.' % lebedev_laikov_npoints)

        if proatomdb is None:
            raise ValueError('When no proatomdb is specified, one can only use the coarse, medium, fine, veryfine, ultrafine or insane grids.')

        atspecs = []
        for n in system.numbers:
            atspecs.append((
                proatomdb.get_rgrid(n),
                nll,
            ))
        return atspecs
