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
'''Rotation of orbitals'''


import numpy as np
from horton.moments import rotate_moments_low, cartesian_transforms


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
    result = np.zeros(coeffs.shape)

    ibasis0 = 0
    for ishell in xrange(obasis.nshell):
        shell_type = obasis.shell_types[ishell]
        if shell_type < 0:
            raise TypeError('Pure functions are not supported in rotate_coeffs.')
        icart0 = ((shell_type+2)*(shell_type+1)*(shell_type))/6
        shellsize = ((shell_type+2)*(shell_type+1))/2
        for ifn in xrange(shellsize):
            rules = cartesian_transforms[icart0+ifn]
            for iorb in xrange(coeffs.shape[1]):
                result[ibasis0+ifn, iorb] = rotate_moments_low(rules, rmat, coeffs[ibasis0:ibasis0+shellsize, iorb])
        ibasis0 += shellsize

    return result
