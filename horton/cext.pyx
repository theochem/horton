# -*- coding: utf-8 -*-
# Horton is a Density Functional Theory program.
# Copyright (C) 2011-2012 Toon Verstraelen <Toon.Verstraelen@UGent.be>
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


import numpy as np
cimport numpy as np
cimport cints
cimport contraction


__all__ = ['fact', 'fact2']


#
# cints wrappers
#

def fact(int n):
    return cints.fact(n)

def fact2(int n):
    return cints.fact2(n)


#
# contraction wrappers
#

def compute_gobasis_overlap(np.ndarray[double, ndim=2] centers,
                            np.ndarray[long, ndim=1] shell_map,
                            np.ndarray[long, ndim=1] num_exponents,
                            np.ndarray[long, ndim=1] num_contractions,
                            np.ndarray[long, ndim=1] con_types,
                            np.ndarray[double, ndim=1] exponents,
                            np.ndarray[double, ndim=1] con_coeffs,
                            np.ndarray[double, ndim=2] output):
    """Compute overlap matrix in Gaussian orbital basis."""
    # TODO: generalize for different storage types of the one-body operator output.
    assert centers.shape[1] == 3
    assert centers.flags['C_CONTIGUOUS']
    assert shell_map.flags['C_CONTIGUOUS']
    assert num_exponents.flags['C_CONTIGUOUS']
    assert num_contractions.flags['C_CONTIGUOUS']
    assert con_types.flags['C_CONTIGUOUS']
    assert exponents.flags['C_CONTIGUOUS']
    assert con_coeffs.flags['C_CONTIGUOUS']
    assert output.shape[0] == output.shape[1]
    assert output.flags['C_CONTIGUOUS']

    # TODO: proper error handling and bounds checking.
    # TODO: provide function pointer for specific integral routine.
    retcode = contraction.compute_gobasis_overlap(
        <double*>centers.data,
        <long*>shell_map.data,
        <long*>num_exponents.data,
        <long*>num_contractions.data,
        <long*>con_types.data,
        <double*>exponents.data,
        <double*>con_coeffs.data,
        <double*>output.data,
    )
    if retcode != 0:
        raise RuntimeError('Something went wrong in compute_gobasis_overlap.')
