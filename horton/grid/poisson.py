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
'''Becke-style numerical Poisson solver'''


import numpy as np

from horton.log import log, timer
from horton.grid.cext import CubicSpline, PotentialExtrapolation
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
            # Derivation of boundary condition at rmax:
            # Multiply differential equation with r**l and integrate. Using
            # partial integration and the fact that V(r)=A/r**(l+1) for large
            # r, we find -(2l+1)A=-4pi*int_0^infty r**2 r**l rho(r) and so
            # V(rmax) = A/rmax**(l+1) = integrate(r**l rho(r))/(2l+1)/rmax**(l+1)
            V_rmax = rgrid.integrate(rho.y*radii**l)/radii[-1]**(l+1)/(2*l+1)
            # Derivation of boundary condition at rmin:
            # Same as for rmax, but multiply differential equation with r**(-l-1)
            # and assume that V(r)=B*r**l for small r.
            V_rmin = rgrid.integrate(rho.y*radii**(-l-1))*radii[0]**(l)/(2*l+1)
            bcs = (V_rmin, None, V_rmax, None)
            v = solve_ode2(b, a, f, bcs, PotentialExtrapolation(l))
            result.append(v)
            counter += 1

    return result
