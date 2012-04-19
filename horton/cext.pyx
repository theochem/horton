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
    'fac2', 'binom', 'gpt_coeff', 'gob_overlap_int1d', 'gob_overlap',
    'get_shell_nbasis', 'gob_normalization',
    'compute_gobasis_overlap', 'I2Gob', 'I2Pow', 'i1pow_inc',
    'project_cartesian_to_pure'
]


#
# cints wrappers
#


def fac2(long n):
    return cints.fac2(n)


def binom(long n, long m):
    return cints.binom(n, m)


def gpt_coeff(long k, long n0, long n1, double pa, double pb):
    return cints.gpt_coeff(k, n0, n1, pa, pb)


def gob_overlap_int1d(long n0, long n1, double pa, double pb, double gamma):
    return cints.gob_overlap_int1d(n0, n1, pa, pb, gamma)


def gob_overlap(double alpha0, long nx0, long ny0, long nz0,
                np.ndarray[double, ndim=1] r0,
                double alpha1, long nx1, long ny1, long nz1,
                np.ndarray[double, ndim=1] r1):
    assert r0.shape[0] == 3
    assert r0.flags['C_CONTIGUOUS']
    assert r1.shape[0] == 3
    assert r1.flags['C_CONTIGUOUS']
    return cints.gob_overlap(alpha0, nx0, ny0, nz0, <double*>r0.data,
                             alpha1, nx1, ny1, nz1, <double*>r1.data)

#
# contraction wrappers
#


def get_shell_nbasis(long shell_type):
    result = contraction.get_shell_nbasis(shell_type)
    if result <= 0:
        raise ValueError("shell_type -1 is not supported.")
    return result


cdef handle_retcode(retcode):
    if retcode == -1:
        raise MemoryError('Ran out of memory inside C routine.')
    elif retcode != 0:
        raise RuntimeError('Unknown error in horton.cext')


cdef check_basis(np.ndarray[double, ndim=2] centers,
                 np.ndarray[long, ndim=1] shell_map,
                 np.ndarray[long, ndim=1] nprims,
                 np.ndarray[long, ndim=1] shell_types,
                 np.ndarray[double, ndim=1] alphas,
                 np.ndarray[double, ndim=1] con_coeffs, long nbasis):
    assert centers.shape[1] == 3
    assert centers.flags['C_CONTIGUOUS']
    assert shell_map.flags['C_CONTIGUOUS']
    assert nprims.flags['C_CONTIGUOUS']
    assert shell_map.shape[0] == nprims.shape[0]
    assert shell_types.flags['C_CONTIGUOUS']
    assert shell_map.shape[0] == shell_types.shape[0]
    assert alphas.flags['C_CONTIGUOUS']
    assert con_coeffs.flags['C_CONTIGUOUS']
    assert alphas.shape[0] == con_coeffs.shape[0]
    assert nprims.sum() == alphas.shape[0]
    assert shell_map.min() >= 0
    assert shell_map.max() < centers.shape[0]
    assert nprims.min() >= 1
    assert (shell_types != -1).any()
    assert nbasis == sum(get_shell_nbasis(shell_type) for shell_type in shell_types)


def compute_gobasis_overlap(np.ndarray[double, ndim=2] centers,
                            np.ndarray[long, ndim=1] shell_map,
                            np.ndarray[long, ndim=1] nprims,
                            np.ndarray[long, ndim=1] shell_types,
                            np.ndarray[double, ndim=1] alphas,
                            np.ndarray[double, ndim=1] con_coeffs,
                            np.ndarray[double, ndim=2] output):
    """Compute overlap matrix in Gaussian orbital basis."""
    # TODO: generalize for different storage types of the one-body operator output.
    assert output.flags['C_CONTIGUOUS']
    assert output.shape[0] == output.shape[1]
    check_basis(centers, shell_map, nprims, shell_types, alphas, con_coeffs, output.shape[0])

    # TODO: provide function pointer for specific integral routine.
    retcode = contraction.compute_gobasis_overlap(
        <double*>centers.data,
        <long*>shell_map.data,
        <long*>nprims.data,
        <long*>shell_types.data,
        <double*>alphas.data,
        <double*>con_coeffs.data,
        shell_map.shape[0],
        centers.shape[0],
        alphas.shape[0],
        <double*>output.data,
        output.shape[0], # nbasis
    )
    handle_retcode(retcode)


def gob_normalization(double alpha, int nx, int ny, int nz):
    return contraction.gob_normalization(alpha, nx, ny, nz)


cdef class I2Gob:
    """Wrapper for the i2gob code, for testing only."""
    cdef contraction.i2gob_type* _c_i2

    # Internal attributes to avoid deallocation.
    cdef np.ndarray _centers
    cdef np.ndarray _shell_map
    cdef np.ndarray _nprims
    cdef np.ndarray _shell_types
    cdef np.ndarray _alphas
    cdef np.ndarray _con_coeffs

    def __cinit__(self):
        self._c_i2 = contraction.i2gob_new()
        if self._c_i2 is NULL:
            raise MemoryError

    def __init__(self, np.ndarray[double, ndim=2] centers,
                 np.ndarray[long, ndim=1] shell_map,
                 np.ndarray[long, ndim=1] nprims,
                 np.ndarray[long, ndim=1] shell_types,
                 np.ndarray[double, ndim=1] alphas,
                 np.ndarray[double, ndim=1] con_coeffs, long nbasis):
        check_basis(centers, shell_map, nprims, shell_types, alphas, con_coeffs, nbasis)

        # Assign arrays to internal attributes
        self._centers = centers
        self._shell_map = shell_map
        self._nprims = nprims
        self._shell_types = shell_types
        self._alphas = alphas
        self._con_coeffs = con_coeffs

        retcode = contraction.i2gob_init(
            self._c_i2,
            <double*>centers.data,
            <long*>shell_map.data,
            <long*>nprims.data,
            <long*>shell_types.data,
            <double*>alphas.data,
            <double*>con_coeffs.data,
            shell_map.shape[0],
            centers.shape[0],
            alphas.shape[0],
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

    def inc_prim(self):
        return contraction.i2gob_inc_prim(self._c_i2)

    def update_prim(self):
        contraction.i2gob_update_prim(self._c_i2)

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

    property alphas:
        def __get__(self):
            return self._alphas

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
                self._c_i2[0].shell_type0, self._c_i2[0].shell_type1,
                self._c_i2[0].alpha0, self._c_i2[0].alpha1,
                self._c_i2[0].r0[0], self._c_i2[0].r0[1], self._c_i2[0].r0[2],
                self._c_i2[0].r1[0], self._c_i2[0].r1[1], self._c_i2[0].r1[2],
                self._c_i2[0].ibasis0, self._c_i2[0].ibasis1,
            )

    property private_fields:
        def __get__(self):
            return (
                self._c_i2[0].ishell0, self._c_i2[0].ishell1,
                self._c_i2[0].nprim0, self._c_i2[0].nprim1,
                self._c_i2[0].oprim0, self._c_i2[0].oprim1,
                self._c_i2[0].iprim0, self._c_i2[0].iprim1,
            )


cdef class I2Pow:
    """Wrapper for the i2pow code, for testing only."""
    cdef contraction.i2pow_type* _c_i2p

    def __cinit__(self):
        self._c_i2p = contraction.i2pow_new()
        if self._c_i2p is NULL:
            raise MemoryError()

    def __init__(self, long shell_type0, long shell_type1, max_nbasis):
        if shell_type0 < 0 or shell_type1 < 0:
            raise ValueError('A shell_type parameter can not be negative.')
        if max_nbasis < get_shell_nbasis(shell_type0):
            raise ValueError('max_nbasis to small for shell_type0.')
        if max_nbasis < get_shell_nbasis(shell_type1):
            raise ValueError('max_nbasis to small for shell_type1.')
        contraction.i2pow_init(self._c_i2p, shell_type0, shell_type1, max_nbasis)

    def inc(self):
        return contraction.i2pow_inc(self._c_i2p)

    def __dealoc__(self):
        if self._c_i2p is not NULL:
            contraction.i2pow_free(self._c_i2p)

    property fields:
        def __get__(self):
            return (
                self._c_i2p[0].nx0, self._c_i2p[0].ny0, self._c_i2p[0].nz0,
                self._c_i2p[0].nx1, self._c_i2p[0].ny1, self._c_i2p[0].nz1,
                self._c_i2p[0].offset
            )


def i1pow_inc(int nx, int ny, int nz):
    result = contraction.i1pow_inc(&nx, &ny, &nz)
    return (nx, ny, nz), result


#
# cartpure wrappers (for testing only)
#

def project_cartesian_to_pure(np.ndarray[double] work_cart,
                              np.ndarray[double] work_pure, long shell_type,
                              long stride, long spacing, long count):
    assert work_cart.flags['C_CONTIGUOUS']
    assert work_pure.flags['C_CONTIGUOUS']
    assert shell_type >= 0
    assert shell_type <= 9
    cartpure.project_cartesian_to_pure(
        <double*> work_cart.data, <double*> work_pure.data, shell_type, stride,
        spacing, count
    )
