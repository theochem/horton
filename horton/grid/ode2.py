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
'''Finite-element second-order ODE solver'''


from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix

from horton.grid.cext import CubicSpline, build_ode2


__all__ = ['solve_ode2']


def solve_ode2(b, a, f, bcs, extrapolation=None):
    '''Solve a second order ODE.

       **Arguments:**

       b, a, f
            Cubic splines for the given functions in the second order ODE. (See
            build_neumann for details.) These cubic splines must have identical
            RTransform objects.

       bcs
            The boundary conditions (See build_neumann for details.)

       **Optional arguments:**

       extrapolation
            The extrapolation object for the returned cubic spline.

       **Returns:** a cubic spline object with the solution that uses the same
       RTransform object as the input functions a, b and f.
    '''
    # Parse args.
    rtf = b.rtransform
    if rtf.to_string() != a.rtransform.to_string():
        raise ValueError('The RTransform objects of b and a do not match.')
    if rtf.to_string() != f.rtransform.to_string():
        raise ValueError('The RTransform objects of b and f do not match.')

    # Transform the given functions to the linear coordinate.
    j1 = rtf.get_deriv()
    j2 = rtf.get_deriv2()
    j3 = rtf.get_deriv3()
    j1sq = j1*j1
    by_new = j1*b.y - j2/j1
    bd_new = j2*b.y + j1sq*b.dx + (j2*j2 - j1*j3)/j1sq
    ay_new = a.y*j1sq
    ad_new = (a.dx*j1sq + 2*a.y*j2)*j1
    fy_new = f.y*j1sq
    fd_new = (f.dx*j1sq + 2*f.y*j2)*j1

    # Transform the boundary conditions
    new_bcs = (
        bcs[0], None if bcs[1] is None else bcs[1]*j1[0],
        bcs[2], None if bcs[3] is None else bcs[3]*j1[-1],
    )

    # Call the equation builder.
    coeffs, rhs = build_ode2(by_new, bd_new, ay_new, ad_new, fy_new, fd_new, new_bcs)
    solution = spsolve(csc_matrix(coeffs), rhs)
    uy_new = solution[::2]
    ud_new = solution[1::2]

    # Transform solution back to the original coordinate.
    uy_orig = uy_new.copy() # A copy of is needed to obtain contiguous arrays.
    ud_orig = ud_new/j1

    return CubicSpline(uy_orig, ud_orig, rtf, extrapolation)
