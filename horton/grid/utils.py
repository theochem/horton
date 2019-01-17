# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2019 The HORTON Development Team
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
'''Auxiliaries for numerical integrals'''


__all__ = ['parse_args_integrate']


def parse_args_integrate(*args, **kwargs):
    '''Parse the arguments given to the integrate method

       **Arguments:**

       data1, data2, ...
            All arguments must be arrays with the same size as the number
            of grid points. The arrays contain the functions, evaluated
            at the grid points, that must be multiplied and integrated.

       **Optional arguments:**

       center=None
            When given, multipole moments are computed with respect to
            this center instead of a plain integral.

       lmax=0
            The maximum angular momentum to consider when computing multipole
            moments

       mtype=1
            The type of multipole moments: 1=``cartesian``, 2=``pure``,
            3=``radial``, 4=``surface``.

       segments=None
            This argument can be used to divide the grid in segments. When
            given, it must be an array with the number of grid points in each
            consecutive segment. The integration is then carried out over each
            segment separately and an array of results is returned. The sum over
            all elements gives back the total integral.

       **Returns:** (args, multipole_args)

       args
            The list of non-keyword arguments, all converted to a flat array.

       multipole_args
            If the center argument is not given, this second return value not
            None. Otherwise, it is a tuple of three elements: center, lmax,
            mtype.

       segments
            The segments argument
    '''
    args = [arg.ravel() for arg in args if arg is not None]

    # process keyword arguments:
    center = kwargs.pop('center', None)
    segments = kwargs.pop('segments', None)
    if center is None:
        # regular integration
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword argument: %s' % kwargs.popitem()[0])

        # Similar to conventional integration routine:
        return args, None, segments
    else:
        # integration with some polynomial factor in the integrand
        lmax = kwargs.pop('lmax', 0)
        mtype = kwargs.pop('mtype', 1)
        if len(kwargs) > 0:
            raise TypeError('Unexpected keyword argument: %s' % kwargs.popitem()[0])
        return args, (center, lmax, mtype), segments
