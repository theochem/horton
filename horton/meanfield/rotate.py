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
'''Rotation of orbitals'''


import numpy as np
from horton.gbasis.cext import fac2
from horton.moments import rotate_cartesian_multipole, get_cartesian_powers


__all__ = ['rotate_coeffs']


def rotate_coeffs(coeffs, obasis, rmat):
    '''Apply a rotation to all cartesian basis functions.

       **Arguments:**

       coeffs
            Expansion coefficients of a set of functions (orbitals) in a local
            basis. shape=(nbasis,nfn)

       obasis
            The local basis set object.

       rmat
            The rotation matrix.
    '''
    if obasis.nbasis != coeffs.shape[0]:
        raise TypeError('The shape of the coefficients array does not match the basis set size')
    if obasis.shell_types.min() < 0:
        raise TypeError('Pure functions are not supported in rotate_coeffs.')

    result = np.zeros(coeffs.shape)

    # 1) undo the part normalization of the basis functions due to the cartesian powers
    lmax = obasis.shell_types.max()
    powers = get_cartesian_powers(lmax)
    factors = []
    for ishell in xrange(obasis.nshell):
        shell_type = obasis.shell_types[ishell]
        icart0 = ((shell_type+2)*(shell_type+1)*(shell_type))/6
        shellsize = ((shell_type+2)*(shell_type+1))/2
        for ifn in xrange(shellsize):
            ipow = icart0+ifn
            factors.append(np.sqrt(fac2(2*powers[ipow,0]-1)*fac2(2*powers[ipow,1]-1)*fac2(2*powers[ipow,2]-1)))
    factors = np.array(factors)
    # replace the array coeffs by the one with undone normalization
    coeffs = coeffs/factors.reshape(-1,1)

    # 2) the actual rotation
    ibasis0 = 0
    for ishell in xrange(obasis.nshell):
        shell_type = obasis.shell_types[ishell]
        icart0 = ((shell_type+2)*(shell_type+1)*(shell_type))/6
        shellsize = ((shell_type+2)*(shell_type+1))/2
        for iorb in xrange(coeffs.shape[1]):
            result[ibasis0:ibasis0+shellsize, iorb] = rotate_cartesian_multipole(rmat, coeffs[ibasis0:ibasis0+shellsize, iorb], 'coeffs')
        ibasis0 += shellsize

    # 3) apply the part of the normalization of the basis functions due to the cartesian powers
    result *= factors.reshape(-1,1)

    return result
