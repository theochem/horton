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
#cython: embedsignature=True
'''C++ extensions for matrix package'''

import numpy as np
cimport numpy as np
np.import_array()

cimport slicing

__all__ = [
    # slicing.cpp
    'compute_slice_abcc','compute_slice_abbc','subtract_slice_abbc',
]

#
# slicing.cpp
#

def compute_slice_abcc(np.ndarray[double, ndim=3] inp not None,
                        np.ndarray[double, ndim=3] inp2 not None,
                        np.ndarray[double, ndim=3] out not None,
                       long nbasis, long nvec):

    assert inp.flags['C_CONTIGUOUS']
    assert inp2.flags['C_CONTIGUOUS']
    assert out.flags['C_CONTIGUOUS']

    assert inp.shape[0] == nvec
    assert inp.shape[1] == inp.shape[2] == nbasis

    assert inp2.shape[0] == nvec
    assert inp2.shape[1] == inp2.shape[2] == nbasis
    assert out.shape[0] == out.shape[1] == out.shape[2] == nbasis

    slicing.get_slice_abcc(&inp[0, 0, 0], &inp2[0, 0, 0], &out[0, 0, 0], nbasis, nvec)
    return out

def compute_slice_abbc(np.ndarray[double, ndim=3] inp not None,
                        np.ndarray[double, ndim=3] inp2 not None,
                        np.ndarray[double, ndim=3] out not None,
                        long nbasis, long nvec):

    assert inp.flags['C_CONTIGUOUS']
    assert inp2.flags['C_CONTIGUOUS']
    assert out.flags['C_CONTIGUOUS']

    assert inp.shape[0] == nvec
    assert inp.shape[1] == inp.shape[2] == nbasis
    assert inp2.shape[0] == nvec
    assert inp2.shape[1] == inp2.shape[2] == nbasis
    assert out.shape[0] == out.shape[1] == out.shape[2] == nbasis

    slicing.get_slice_abbc(&inp[0, 0, 0], &inp2[0, 0, 0], &out[0, 0, 0], nbasis,nvec)
    return out

def subtract_slice_abbc(np.ndarray[double, ndim=3] inp not None,
                         np.ndarray[double, ndim=3] inp2 not None,
                         np.ndarray[double, ndim=3] out not None,
                         long nbasis, long nvec):

    assert inp.flags['C_CONTIGUOUS']
    assert inp2.flags['C_CONTIGUOUS']
    assert out.flags['C_CONTIGUOUS']

    assert inp.shape[0] == nvec
    assert inp.shape[1] == inp.shape[2] == nbasis
    assert inp2.shape[0] == nvec
    assert inp2.shape[1] == inp2.shape[2] == nbasis
    assert out.shape[0] == out.shape[1] == out.shape[2] == nbasis

    slicing.sub_slice_abbc(&inp[0, 0, 0], &inp2[0, 0, 0], &out[0, 0, 0], nbasis,nvec)
    return out
