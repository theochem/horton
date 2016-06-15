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
'''C++ extensions for matrix package'''

import numpy as np
cimport numpy as np
np.import_array()

cimport libcpp

cimport slicing

__all__ = [
    # slicing.cpp
    'slice_to_three_abbc_abc', 'slice_to_three_abcc_bac',
    'slice_to_three_abcc_abc',
]

#
# slicing.cpp
#

def _check_args_slice_to_three(np.ndarray[double, ndim=3] inp not None,
                               np.ndarray[double, ndim=3] inp2 not None,
                               np.ndarray[double, ndim=3] out not None):
    '''Check the arguments of the functions below'''
    assert inp.flags['C_CONTIGUOUS']
    assert inp2.flags['C_CONTIGUOUS']
    assert out.flags['C_CONTIGUOUS']
    cdef long nvec = inp.shape[0]
    cdef long nbasis = inp.shape[1]
    assert inp.shape[2] == nbasis
    assert inp2.shape[0] == nvec
    assert inp2.shape[1] == inp2.shape[2] == nbasis
    assert out.shape[0] == out.shape[1] == out.shape[2] == nbasis


def slice_to_three_abbc_abc(np.ndarray[double, ndim=3] inp not None,
                            np.ndarray[double, ndim=3] inp2 not None,
                            np.ndarray[double, ndim=3] out not None,
                            double factor=1.0, libcpp.bool clear=True):
    _check_args_slice_to_three(inp, inp2, out)
    slicing.slice_to_three_abbc_abc(&inp[0, 0, 0], &inp2[0, 0, 0], &out[0, 0, 0], factor, clear, inp.shape[1], inp.shape[0])



def slice_to_three_abcc_bac(np.ndarray[double, ndim=3] inp not None,
                            np.ndarray[double, ndim=3] inp2 not None,
                            np.ndarray[double, ndim=3] out not None,
                            double factor=1.0, libcpp.bool clear=True):
    _check_args_slice_to_three(inp, inp2, out)
    slicing.slice_to_three_abcc_bac(&inp[0, 0, 0], &inp2[0, 0, 0], &out[0, 0, 0], factor, clear, inp.shape[1], inp.shape[0])


def slice_to_three_abcc_abc(np.ndarray[double, ndim=3] inp not None,
                            np.ndarray[double, ndim=3] inp2 not None,
                            np.ndarray[double, ndim=3] out not None,
                            double factor=1.0, libcpp.bool clear=True):
    _check_args_slice_to_three(inp, inp2, out)
    slicing.slice_to_three_abcc_abc(&inp[0, 0, 0], &inp2[0, 0, 0], &out[0, 0, 0], factor, clear, inp.shape[1], inp.shape[0])
