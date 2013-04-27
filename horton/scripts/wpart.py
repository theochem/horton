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


from horton.grid.int1d import SimpsonIntegrator1D
from horton.grid.cext import lebedev_laikov_npoints

__all__ = ['parse_grid']


def parse_grid(sgrid, system, proatomdb):
    if sgrid == 'coarse':
        return 'tv-13.1-3'
    elif sgrid == 'medium':
        return 'tv-13.1-4'
    elif sgrid == 'fine':
        return 'tv-13.1-5'
    elif sgrid == 'veryfine':
        return 'tv-13.1-6'
    else:
        try:
            nll = int(sgrid)
        except ValueError:
            raise ValueError('If the grid argument is not coarse, medium, fine or veryfine, it should be an integer.')

        if nll not in lebedev_laikov_npoints:
            raise ValueError('If the grid argument is an integer, it should be one of the following: %s.' % lebedev_laikov_npoints)

        atspecs = []
        for n in system.numbers:
            atspecs.append((
                proatomdb.get_rgrid(n).rtransform,
                SimpsonIntegrator1D(),
                nll,
            ))
        return atspecs
