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

cimport cartpure
cimport cints
cimport contraction


__all__ = [
    'fact', 'fact2', 'cints_overlap',
    'get_con_nbasis', 'gob_normalization',
    'compute_gobasis_overlap', 'I2Gob', 'I2Pow', 'i1pow_inc',
    'project_cartesian_to_pure'
]


#
# cints wrappers
#

def fact(int n):
    return cints.fact(n)

def fact2(int n):
    return cints.fact2(n)

def cints_overlap(double alpha1, int l1, int m1, int n1,
                  double xa, double ya, double za,
                  double alpha2, int l2, int m2, int n2,
                  double xb, double yb, double zb):
    return cints.overlap(alpha1, l1, m1, n1, xa, ya, za, alpha2, l2, m2, n2, xb,
                         yb, zb)

#
# contraction wrappers
#


def get_con_nbasis(long con_type):
    result = contraction.get_con_nbasis(con_type)
    if result <= 0:
        raise ValueError("con_type -1 is not supported.")
    return result


cdef handle_retcode(retcode):
    if retcode == -1:
        raise MemoryError('Ran out of memory inside C routine.')
    elif retcode == -2:
        raise ValueError('The center indexes in the shell_map are out of range.')
    elif retcode == -3:
        raise ValueError('The elements of ncons should be at least 1.')
    elif retcode == -4:
        raise ValueError('The size of the array con_types does not match the sum of ncons.')
    elif retcode == -5:
        raise ValueError('The elements of nexps should be at least 1.')
    elif retcode == -6:
        raise ValueError('The size of the array exponents does not match the sum of nexps.')
    elif retcode == -7:
        raise ValueError('Encountered the nonexistent con_type -1.')
    elif retcode == -8:
        raise ValueError('The size of con_coeffs does not match nexp and ncon.')
    elif retcode == -9:
        raise ValueError('The total number of basis functions does not match the output array.')
    elif retcode < 0:
        raise RuntimeError('Unknown error in horton.cext')


def compute_gobasis_overlap(np.ndarray[double, ndim=2] centers,
                            np.ndarray[long, ndim=1] shell_map,
                            np.ndarray[long, ndim=1] ncons,
                            np.ndarray[long, ndim=1] nexps,
                            np.ndarray[long, ndim=1] con_types,
                            np.ndarray[double, ndim=1] exponents,
                            np.ndarray[double, ndim=1] con_coeffs,
                            np.ndarray[double, ndim=2] output):
    """Compute overlap matrix in Gaussian orbital basis."""
    # TODO: generalize for different storage types of the one-body operator output.
    assert centers.shape[1] == 3
    assert centers.flags['C_CONTIGUOUS']
    assert shell_map.flags['C_CONTIGUOUS']
    assert ncons.flags['C_CONTIGUOUS']
    assert nexps.shape[0] == ncons.shape[0]
    assert nexps.flags['C_CONTIGUOUS']
    assert shell_map.shape[0] == nexps.shape[0]
    assert con_types.flags['C_CONTIGUOUS']
    assert exponents.flags['C_CONTIGUOUS']
    assert con_coeffs.flags['C_CONTIGUOUS']
    assert output.shape[0] == output.shape[1]
    assert output.flags['C_CONTIGUOUS']

    # TODO: provide function pointer for specific integral routine.
    retcode = contraction.compute_gobasis_overlap(
        <double*>centers.data,
        <long*>shell_map.data,
        <long*>ncons.data,
        <long*>nexps.data,
        <long*>con_types.data,
        <double*>exponents.data,
        <double*>con_coeffs.data,
        len(shell_map),
        len(centers),
        len(con_types),
        len(exponents),
        len(con_coeffs),
        <double*>output.data,
        len(output), # nbasis
    )
    handle_retcode(retcode)


def gob_normalization(double alpha, int l, int m, int n):
    return contraction.gob_normalization(alpha, l, m, n)


cdef class I2Gob:
    """Wrapper for the i2gob code, for testing only."""
    cdef contraction.i2gob_type* _c_i2

    # Internal attributes to avoid deallocation.
    cdef np.ndarray _centers
    cdef np.ndarray _shell_map
    cdef np.ndarray _ncons
    cdef np.ndarray _nexps
    cdef np.ndarray _con_types
    cdef np.ndarray _exponents
    cdef np.ndarray _con_coeffs

    def __cinit__(self):
        self._c_i2 = contraction.i2gob_new()
        if self._c_i2 is NULL:
            raise MemoryError

    def __init__(self, np.ndarray[double, ndim=2] centers,
                 np.ndarray[long, ndim=1] shell_map,
                 np.ndarray[long, ndim=1] ncons,
                 np.ndarray[long, ndim=1] nexps,
                 np.ndarray[long, ndim=1] con_types,
                 np.ndarray[double, ndim=1] exponents,
                 np.ndarray[double, ndim=1] con_coeffs, long nbasis):

        assert centers.shape[1] == 3
        assert centers.flags['C_CONTIGUOUS']
        assert shell_map.flags['C_CONTIGUOUS']
        assert ncons.flags['C_CONTIGUOUS']
        assert nexps.shape[0] == ncons.shape[0]
        assert nexps.flags['C_CONTIGUOUS']
        assert shell_map.shape[0] == nexps.shape[0]
        assert con_types.flags['C_CONTIGUOUS']
        assert exponents.flags['C_CONTIGUOUS']
        assert con_coeffs.flags['C_CONTIGUOUS']
        assert nbasis > 0

        # Assign arrays to internal attributes
        self._centers = centers
        self._shell_map = shell_map
        self._ncons = ncons
        self._nexps = nexps
        self._con_types = con_types
        self._exponents = exponents
        self._con_coeffs = con_coeffs

        # TODO: proper error handling and bounds checking.
        # TODO: provide function pointer for specific integral routine.
        retcode = contraction.i2gob_init(
            self._c_i2,
            <double*>centers.data,
            <long*>shell_map.data,
            <long*>ncons.data,
            <long*>nexps.data,
            <long*>con_types.data,
            <double*>exponents.data,
            <double*>con_coeffs.data,
            len(shell_map),
            len(centers),
            len(con_types),
            len(exponents),
            len(con_coeffs),
            nbasis,
        )
        handle_retcode(retcode)

    def __dealoc__(self):
        if self._c_i2 is not NULL:
            contraction.i2gob_free(self._c_i2)

    def inc_shell(self):
        return contraction.i2gob_inc_shell(self._c_i2)

    def update_shell(self):
        contraction.i2gob_update_shell(self._c_i2)

    def inc_con(self):
        return contraction.i2gob_inc_con(self._c_i2)

    def update_con(self):
        contraction.i2gob_update_con(self._c_i2)

    def inc_exp(self):
        return contraction.i2gob_inc_exp(self._c_i2)

    def update_exp(self):
        contraction.i2gob_update_exp(self._c_i2)

    def store(self, np.ndarray[double, ndim=2] work_pure,
              np.ndarray[double, ndim=2] output):
        assert work_pure.shape[0] == self._c_i2[0].max_nbasis
        assert work_pure.shape[1] == self._c_i2[0].max_nbasis
        assert work_pure.flags['C_CONTIGUOUS']
        assert output.shape[0] == self._c_i2[0].nbasis
        assert output.shape[1] == self._c_i2[0].nbasis
        assert output.flags['C_CONTIGUOUS']
        contraction.i2gob_store(self._c_i2, <double*>work_pure.data, <double*>output.data)

    property centers:
        def __get__(self):
            return self._centers

    property exponents:
        def __get__(self):
            return self._exponents

    property con_coeffs:
        def __get__(self):
            return self._con_coeffs

    property max_nbasis:
        def __get__(self):
            return self._c_i2[0].max_nbasis

    property public_fields:
        def __get__(self):
            return (
                self._c_i2[0].con_coeff,
                self._c_i2[0].con_type0, self._c_i2[0].con_type1,
                self._c_i2[0].exp0, self._c_i2[0].exp1,
                self._c_i2[0].x0, self._c_i2[0].y0, self._c_i2[0].z0,
                self._c_i2[0].x1, self._c_i2[0].y1, self._c_i2[0].z1,
                self._c_i2[0].ibasis0, self._c_i2[0].ibasis1,
            )

    property private_fields:
        def __get__(self):
            return (
                self._c_i2[0].ishell0, self._c_i2[0].ishell1,
                self._c_i2[0].ncon0, self._c_i2[0].ncon1,
                self._c_i2[0].ocon0, self._c_i2[0].ocon1,
                self._c_i2[0].icon0, self._c_i2[0].icon1,
                self._c_i2[0].nexp0, self._c_i2[0].nexp1,
                self._c_i2[0].oexp0, self._c_i2[0].oexp1,
                self._c_i2[0].iexp0, self._c_i2[0].iexp1,
                self._c_i2[0].occ0, self._c_i2[0].occ1,
            )


cdef class I2Pow:
    """Wrapper for the i2pow code, for testing only."""
    cdef contraction.i2pow_type* _c_i2p

    def __cinit__(self):
        self._c_i2p = contraction.i2pow_new()
        if self._c_i2p is NULL:
            raise MemoryError()

    def __init__(self, long con_type0, long con_type1, max_nbasis):
        if con_type0 < 0 or con_type1 < 0:
            raise ValueError('A con_type parameter can not be negative.')
        if max_nbasis < get_con_nbasis(con_type0):
            raise ValueError('max_nbasis to small for con_type0.')
        if max_nbasis < get_con_nbasis(con_type1):
            raise ValueError('max_nbasis to small for con_type1.')
        contraction.i2pow_init(self._c_i2p, con_type0, con_type1, max_nbasis)

    def inc(self):
        return contraction.i2pow_inc(self._c_i2p)

    def __dealoc__(self):
        if self._c_i2p is not NULL:
            contraction.i2pow_free(self._c_i2p)

    property fields:
        def __get__(self):
            return (
                self._c_i2p[0].l0, self._c_i2p[0].m0, self._c_i2p[0].n0,
                self._c_i2p[0].l1, self._c_i2p[0].m1, self._c_i2p[0].n1,
                self._c_i2p[0].offset
            )


def i1pow_inc(int l, int m, int n):
    result = contraction.i1pow_inc(&l, &m, &n)
    return (l, m, n), result


#
# cartpure wrappers (for testing only)
#

def project_cartesian_to_pure(np.ndarray[double] work_cart,
                              np.ndarray[double] work_pure, long con_type,
                              long stride, long spacing, long count):
    assert work_cart.flags['C_CONTIGUOUS']
    assert work_pure.flags['C_CONTIGUOUS']
    assert con_type >= 0
    assert con_type <= 9
    cartpure.project_cartesian_to_pure(
        <double*> work_cart.data, <double*> work_pure.data, con_type, stride,
        spacing, count
    )
