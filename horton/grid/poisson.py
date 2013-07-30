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
'''Becke-style numerical Poisson solver'''


import numpy as np

from horton.log import log, timer
from horton.grid.cext import CubicSpline, PowerExtrapolation
from horton.grid.ode2 import solve_ode2
from horton.grid.radial import RadialGrid


__all__ = ['solve_poisson_becke']


@timer.with_section('Becke Poisson')
def solve_poisson_becke(density_decomposition):
    '''Compute the electrostatic potential of a density expanded in real spherical harmonics

       **Arguments:**

       density_decomposition
            A list of cubic splines returned by the method
            AtomicGrid.get_spherical_decomposition.

       The returned list of splines is a spherical decomposition of the
       hartree potential (felt by a particle with the same charge unit as the
       density).
    '''
    log.cite('becke1988_poisson', 'the numerical integration of the Poisson equation')

    lmax = np.sqrt(len(density_decomposition)) - 1
    assert lmax == int(lmax)
    lmax = int(lmax)

    result = []
    counter = 0
    for l in xrange(0, lmax+1):
        for m in xrange(-l, l+1):
            rho = density_decomposition[counter]
            rtf = rho.rtransform
            rgrid = RadialGrid(rtf)
            radii = rtf.get_radii()
            # The approach followed here is obtained after substitution of
            # u = r*V in Eq. (21) in Becke's paper. After this transformation,
            # the boundary conditions can be implemented such that the output
            # is more accurate.
            fy = -4*np.pi*rho.y
            fd = -4*np.pi*rho.dx
            f = CubicSpline(fy, fd, rtf)
            b = CubicSpline(2/radii, -2/radii**2, rtf)
            a = CubicSpline(-l*(l+1)*radii**-2, 2*l*(l+1)*radii**-3, rtf)
            bcs = (None, 0.0, rgrid.integrate(rho.y)/radii[-1]**(l+1), None)
            v = solve_ode2(b, a, f, bcs, PowerExtrapolation(-l-1))
            result.append(v)
            counter += 1

    return result
